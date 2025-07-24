function record = analyse_njtestrecord(record,verbose)
%analyse_njtestrecord Analyses Nora Jamann's Neuropixels experiment
%
%  RECORD = analyse_njtestrecord( RECORD, [VERBOSE=true])
%
% 2024, Alexander Heimel

if nargin<2 || isempty(verbose)
    verbose = true;
end

params = nj_default_parameters();

if params.use_fixed_random_seed
    rng(19750321);
end

sAP = nj_load_data(record);

[spont_start,spont_end,spont_duration] = nj_get_spont_interval(record,sAP);
pulse_times = nj_get_pulse_times(record,sAP);
n_pulses = sum(cellfun(@length,pulse_times));

if params.use_only_good_clusters
    % very stringent:
    % good_clusters = strcmp({sAP.sCluster.bc_unitType},'GOOD'); % only somatic good

    % stringent:
    % good_clusters = contains({sAP.sCluster.bc_unitType}, 'GOOD'); % includes non-somatic

    % less stringent:
    good_clusters = contains({sAP.sCluster.bc_unitType}, 'GOOD') | ([sAP.sCluster.KilosortGood]==1);
 
    % note: selection of only GOOD units is done later in
    % make_figures_jamann_neuropixels

    sAP.sCluster = sAP.sCluster(good_clusters);
end

% select only barrel cortex layer 5 and thalamus
ind = (contains({sAP.sCluster.Area},'barrel') & contains({sAP.sCluster.Area},'layer 5')) | contains({sAP.sCluster.Area},'thalamus');

% select cortex and thalamus
%ind = (contains({sAP.sCluster.Area},'layer')) | contains({sAP.sCluster.Area},'thalamus');

sAP.sCluster = sAP.sCluster(ind);
n_clusters = length(sAP.sCluster);

%
%sAP
%   sCluster - clusters
%      SpikeTimes - [n x 1] all spiketims
%      Waveforn - [1 x n] Waveform
%      Area - char array with area name
%      KilosortGood - 0 or 1
%      KilosortLabel - 'mua', ...
%      bc_unitType - BombCell type 'MUA', 'GOOD', ...

% Compute spike crosscorrelation
measures = [];

tbin_centers = -0.100:params.dt:0.100;
xcortbins = [min(tbin_centers)-params.dt/2 tbin_centers+params.dt/2];

period_types = {'all','spont_test','opto_test','whisker_test','opto_plus_whisker'};


n_period_types= length(period_types);

h_wait = waitbar(0,'Please wait...');

m = 0;


for i = 1:n_clusters

    if params.remove_spike_doublets
        d = diff(sAP.sCluster(i).SpikeTimes);
        ind = find(d==0,1);
        if ~isempty(ind)
            ind = find(d==0);
            logmsg(['Removed ' num2str(length(ind)) ' spike doublets in cluster ' num2str(i)]);
            sAP.sCluster(i).SpikeTimes(ind) = [];
        end
    end

    for p = 1:n_period_types
        period_type = period_types{p};

        m = m + 1;
        measures(m).cluster = sAP.sCluster(i).Cluster;
        measures(m).area = sAP.sCluster(i).Area;
        measures(m).good_bc = contains(sAP.sCluster(i).bc_unitType,'GOOD'); % include non-somatic good
        measures(m).good_ks = sAP.sCluster(i).KilosortGood;
        measures(m).non_somatic = contains(sAP.sCluster(i).bc_unitType,'NON-SOMA');
        measures(m).period_type = period_type;


        switch period_type
            case 'all'
                sp_i = sAP.sCluster(i).SpikeTimes;
                duration = sp_i(end)-sp_i(1);
            case {'opto_test','whisker_test','opto_plus_whisker'}

                xcorrv = zeros(length(tbin_centers),1);
                duration = 0;
                ind_all = [];
                n_pulses = 0;
                for t = record.(period_type)(:)'
                    if isempty(pulse_times{t})
                        continue
                    end
                    start_time = pulse_times{t}(1) - 0.1;
                    end_time = pulse_times{t}(end) + 0.1;
                    duration = duration + end_time - start_time;
                    n_pulses = n_pulses + length(pulse_times{t});
                    ind = find(sAP.sCluster(i).SpikeTimes>start_time & sAP.sCluster(i).SpikeTimes<end_time) ;
                    if isempty(ind)
                        continue
                    end
                    ind_all = [ind_all;ind];
                    xcorrv = xcorrv + xcorrspiketimes(pulse_times{t},sAP.sCluster(i).SpikeTimes(ind),xcortbins);
                end
                sp_i = sAP.sCluster(i).SpikeTimes(ind_all);

                ind_post = (tbin_centers>0 & tbin_centers<params.max_delay);
                n_post_spikes =  sum(xcorrv(ind_post));
                measures(m).pulse_transfer = n_post_spikes / n_pulses ;
                measures(m).zeta_p = min([sAP.sCluster(i).ZetaP(record.(period_type))]);
            case 'spont_test'
                ind = find(sAP.sCluster(i).SpikeTimes>spont_start & sAP.sCluster(i).SpikeTimes<spont_end) ;
                sp_i = sAP.sCluster(i).SpikeTimes(ind);
        end

        [bursts_i,tonics_i,n_events_i,spikes_per_burst_i,mean_burst_isi_i] = get_bursttimes(sp_i,params);

        measures(m).rate = length(sp_i)/duration;
        measures(m).burst_rate = length(bursts_i)/duration;
        measures(m).n_spikes_per_burst = mean(spikes_per_burst_i);
        measures(m).burst_fraction = length(bursts_i) / n_events_i;

        if params.bursting_threshold<1
            measures(m).bursting = measures(i).burst_fraction > params.bursting_threshold;
        else
            measures(m).bursting = (length(bursts_i)>=params.bursting_threshold); % Nora's definition
        end

        measures(m).burst_frequency = 1/mean_burst_isi_i;


        if length(sp_i)>1
            measures(m).isi_hist_edges = params.isi_hist_edges;
            measures(m).isi_hist_frac = histcounts(diff(sp_i),'BinEdges',params.isi_hist_edges)/(length(sp_i)-1);
            mask_below_1s = (params.isi_hist_edges<=1);
            measures(m).frac_below_1s = sum(measures(m).isi_hist_frac(mask_below_1s));
       else
            measures(m).isi_hist_edges = [];
            measures(m).isi_hist_frac = [];
            measures(m).frac_below_1s = NaN;
        end
 
        measures(m).pairs = [];
        if isempty(sp_i) 
            continue
        end
        if ~params.analyze_pairs
            continue
        end
        if ~contains(sAP.sCluster(i).Area,'layer' ) % only look at cortex as presynaptic neuron
            continue
        end
        pairs = get_pairs( sAP, period_type, pulse_times, sp_i, bursts_i, tonics_i, spont_start,spont_end,tbin_centers,xcortbins,record,params);
        measures(m).pairs = pairs;
    end % period_type p
    waitbar(i/n_clusters,h_wait);
end % cluster i
close(h_wait)

record.measures = measures;

% pick best period and add as measures
record = nj_combine_periods(record);

logmsg('Done')

end




function pairs = get_pairs( sAP, period_type, pulse_times, sp_i, bursts_i, tonics_i, spont_start,spont_end,tbin_centers,xcortbins,record,params  )
% find cells that are postsynaptically connected based on spike cross correlogram

    

% look for connected pairs
count = 0;
pairs = struct([]);

if length(sp_i)<4 % too few spikes to say anything
    return
end



first_spikes_i = union(bursts_i,tonics_i);

n_clusters = length(sAP.sCluster);


% only look at first spikes in burst
if params.only_pair_on_first_spikes
    if isempty(first_spikes_i)
        keyboard
    end

    
    % if true, then burst transfer should still be calculated differently
    sp_i = first_spikes_i;

end

for j = 1:n_clusters

    if ~contains(sAP.sCluster(j).Area,'thalamus' )
        continue
    end

    switch period_type
        case 'all'
            sp_j = sAP.sCluster(j).SpikeTimes;
        case {'opto_test','whisker_test','opto_plus_whisker'}
            ind_all = [];
            for t = record.(period_type)(:)'
                if isempty(pulse_times{t})
                    continue
                end
                start_time = pulse_times{t}(1) - 0.1;
                end_time = pulse_times{t}(end) + 0.1;
                ind = find(sAP.sCluster(j).SpikeTimes>start_time & sAP.sCluster(j).SpikeTimes<end_time) ;
                ind_all = [ind_all;ind];
            end
            sp_j = sAP.sCluster(j).SpikeTimes(ind_all);
        case 'spont_test'
            ind = find(sAP.sCluster(j).SpikeTimes>spont_start & sAP.sCluster(j).SpikeTimes<spont_end) ;
            sp_j = sAP.sCluster(j).SpikeTimes(ind);
    end

    if isempty(sp_j)
        continue
    end

    xcorrv = xcorrspiketimes(sp_i,sp_j,xcortbins);

    [m,ind_peak] = max(xcorrv);
    peak_time = tbin_centers(ind_peak);

    if peak_time<=0 % only get pre before post
        continue
    end
    if abs(peak_time) > params.max_delay   % peak has to be within 20 ms
        continue
    end
    if m < mean(xcorrv) + 3 * std(xcorrv) % peak has to be larger than average
        continue
    end
    if m < 4 % less than 4 spikes in total will be too little to be significant
        continue
    end

    xcorrv_jit = NaN([params.n_jitters size(tbin_centers)]);
    for k = 1:params.n_jitters
        sp1_jit = sp_i + 0.02 *(rand(size(sp_i))-0.5);
        sp2_jit = sp_j + 0.02 *(rand(size(sp_j))-0.5);
        xcorrv_jit(k,:) = xcorrspiketimes(sp1_jit,sp2_jit,xcortbins);
    end %k
    xcorrv_mean = squeeze(mean(xcorrv_jit,1));
    xcorrv_std = squeeze(std(xcorrv_jit,1));

    peak_snr = (xcorrv(ind_peak) - max(xcorrv_mean)) / max(xcorrv_std);

    if peak_snr < 3 % to limit analyzing too many peaks
        continue
    end

    % if we are here, we have discovered a putative connection
    count = count + 1;

    pairs(count).peak_snr = peak_snr;
    pairs(count).peak_time = peak_time;
    pairs(count).peak_height = xcorrv(ind_peak);

    % fraction of first 10 ms post over first 20 ms post
    ind_5ms = (tbin_centers>0 & tbin_centers<=0.005);
    ind_10ms = (tbin_centers>0 & tbin_centers<=0.010);
    ind_20ms = (tbin_centers>0 & tbin_centers<=0.020);
    if sum(xcorrv(ind_20ms))>0
        pairs(count).fraction_post_10ms = sum(xcorrv(ind_10ms))/sum(xcorrv(ind_20ms));
        pairs(count).fraction_post_5ms = sum(xcorrv(ind_5ms))/sum(xcorrv(ind_20ms));
    else
        pairs(count).fraction_post_10ms = NaN;
         pairs(count).fraction_post_5ms = NaN;
   end
 
    % time of first bin post presynaptic spike with non-zero spikes
    ind_first_nonzero_bin = find( xcorrv(:)>0 & tbin_centers(:)>0,1);
    pairs(count).delay_first_spikes = tbin_centers(ind_first_nonzero_bin);

    % fraction of first 10 ms post first_nonzero_bin over 20 ms post
    ind_5ms = (tbin_centers>=pairs(count).delay_first_spikes & tbin_centers<=pairs(count).delay_first_spikes+0.005);
    ind_10ms = (tbin_centers>=pairs(count).delay_first_spikes & tbin_centers<=pairs(count).delay_first_spikes+0.010);
    ind_20ms = (tbin_centers>=pairs(count).delay_first_spikes & tbin_centers<=pairs(count).delay_first_spikes+0.020);
    if sum(xcorrv(ind_20ms))>0
        pairs(count).fraction_10ms_post_first = sum(xcorrv(ind_10ms))/sum(xcorrv(ind_20ms));
        pairs(count).fraction_5ms_post_first = sum(xcorrv(ind_5ms))/sum(xcorrv(ind_20ms));
    else
        pairs(count).fraction_5ms_post_first = NaN;
        pairs(count).fraction_10ms_post_first = NaN;
    end

    % time of first bin post presynaptic spike with more spikes than before presynaptic spike
    peak_before = max( xcorrv(tbin_centers<0));
    ind_first_nonzero_bin = find( xcorrv(:)>peak_before & tbin_centers(:)>0,1);
    pairs(count).delay_first_increased_spikes = tbin_centers(ind_first_nonzero_bin);

    % fraction of first 10 ms post increased bin over 20 ms post
    ind_5ms = (tbin_centers>=pairs(count).delay_first_increased_spikes & tbin_centers<=pairs(count).delay_first_increased_spikes+0.005);
    ind_10ms = (tbin_centers>=pairs(count).delay_first_increased_spikes & tbin_centers<=pairs(count).delay_first_increased_spikes+0.010);
    ind_20ms = (tbin_centers>=pairs(count).delay_first_increased_spikes & tbin_centers<=pairs(count).delay_first_increased_spikes+0.020);
    if sum(xcorrv(ind_20ms))>0
        pairs(count).fraction_10ms_post_increase = sum(xcorrv(ind_10ms))/sum(xcorrv(ind_20ms));
        pairs(count).fraction_5ms_post_increase = sum(xcorrv(ind_5ms))/sum(xcorrv(ind_20ms));
    else
        pairs(count).fraction_5ms_post_increase = NaN;
        pairs(count).fraction_10ms_post_increase = NaN;
    end



    [bursts_j,tonics_j] = get_bursttimes(sp_j,params);
    first_spikes_j = union(bursts_j,tonics_j);

    xcorrv_first_spikes = xcorrspiketimes(first_spikes_i,first_spikes_j,xcortbins);

    if ~isempty(bursts_i) && ~isempty(bursts_j)
        xcorrv_bursts = xcorrspiketimes(bursts_i,bursts_j,xcortbins);
    else
        xcorrv_bursts = zeros(size(tbin_centers));
    end

    if ~isempty(bursts_i) && ~isempty(sp_j)
        xcorrv_bursts_spikes = xcorrspiketimes(bursts_i,sp_j,xcortbins);
    else
        xcorrv_bursts_spikes =zeros(size(tbin_centers));
    end

    ind_post = (tbin_centers>0 & tbin_centers<params.max_delay);
    n_post_bursts =  sum(xcorrv_bursts(ind_post));
    pairs(count).burst_transfer = n_post_bursts / length(bursts_i);

    n_post_spikes = sum(xcorrv(ind_post));
    pairs(count).spike_transfer = n_post_spikes / length(sp_i);

    n_spikes_post_burst = sum(xcorrv_bursts_spikes(ind_post));
    pairs(count).spikes_post_burst = n_spikes_post_burst / length(bursts_i);

    
    if params.compute_individual_burst_success
        n_bursts = length(bursts_i);
        burst_isi = NaN(1,n_bursts);
        burst_transferred = NaN(1,n_bursts);
        for burst = 1:n_bursts
            burst_time = bursts_i(burst);

            ind_spikes_after_burst = find(sp_j>burst_time & sp_j<burst_time+params.max_delay,1);
            burst_transferred(burst) = (~isempty(ind_spikes_after_burst));

            ind_spike = find(sp_i==burst_time);
            if length(ind_spike)>1
                disp(['Spike doublet at t = ' num2str(burst_time)])
                ind_spike = ind_spike(1);
            end
            burst_isi(burst) = sp_i(ind_spike+1)-sp_i(ind_spike);
        end
        if any(burst_transferred)
            pairs(count).transferred_burst_isi = mean(burst_isi(logical(burst_transferred)),'omitnan');
        else
            pairs(count).transferred_burst_isi = NaN;
        end
        if ~all(burst_transferred)
            pairs(count).not_transferred_burst_isi = mean(burst_isi(~logical(burst_transferred)),'omitnan');
        else
            pairs(count).not_transferred_burst_isi = NaN;
        end
        fast_bursts_mask = burst_isi<params.fast_burst_max_isi;
        slow_bursts_mask = ~fast_bursts_mask;

        pairs(count).fast_burst_transfer = sum(fast_bursts_mask & burst_transferred)/sum(fast_bursts_mask);
        pairs(count).slow_burst_transfer = sum(slow_bursts_mask & burst_transferred)/sum(slow_bursts_mask);
        
    end


    pairs(count).post_spike_time_mean = (tbin_centers(ind_post)*xcorrv(ind_post))/sum(xcorrv(ind_post));
    v = ((tbin_centers(ind_post).^2)*xcorrv(ind_post))/sum(xcorrv(ind_post)) - ...
        ((tbin_centers(ind_post)*xcorrv(ind_post))/sum(xcorrv(ind_post)))^2 ;
    pairs(count).post_spike_time_std = sqrt(v);

    % first spike mean, variance and std 
    pairs(count).post_first_spike_time_mean = (tbin_centers(ind_post)*xcorrv_first_spikes(ind_post))/sum(xcorrv_first_spikes(ind_post));
    v = ((tbin_centers(ind_post).^2)*xcorrv_first_spikes(ind_post))/sum(xcorrv_first_spikes(ind_post)) - ...
        ((tbin_centers(ind_post)*xcorrv_first_spikes(ind_post))/sum(xcorrv_first_spikes(ind_post)))^2 ;
    pairs(count).post_first_spike_time_var = v;
    pairs(count).post_first_spike_time_std = sqrt(v);


    % pairs(count).pre_cluster = sAP.sCluster(i).Cluster; %#ok<*AGROW>
    % pairs(count).pre_area = sAP.sCluster(i).Area;
    pairs(count).post_cluster = sAP.sCluster(j).Cluster;
    pairs(count).post_area = sAP.sCluster(j).Area;
    pairs(count).tbin_centers = tbin_centers;
    pairs(count).xcorrv = xcorrv;
    pairs(count).xcorrv_mean = xcorrv_mean;
    pairs(count).xcorrv_std = xcorrv_std;
    pairs(count).xcorrv_first_spikes = xcorrv_first_spikes;
    pairs(count).xcorrv_bursts = xcorrv_bursts;

    pairs(count).post_cluster_good_bc = strcmp(sAP.sCluster(j).bc_unitType,'GOOD');
    pairs(count).post_cluster_good_ks = sAP.sCluster(j).KilosortGood;
end % j



end


