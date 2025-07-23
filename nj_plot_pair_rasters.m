function nj_plot_pair_rasters(record,pre_cluster,post_cluster,period_type,sAP)
%nj_plot_pair_rasters. Plot raster of pre- and post-synaptic pair
% 
%   nj_plot_pair_rasters(record,pre_cluster,post_cluster,period_type,sAP)
%
% 2025, Alexander Heimel


if nargin<5 || isempty(sAP)
    sAP = nj_load_data(record);
end
if nargin<4 || isempty(period_type)
    period_type = '';
end

params = nj_default_parameters(record);

if params.use_fixed_random_seed
    rng(1);
end

measures = record.measures;
sel = ['cluster=' num2str(pre_cluster)];
if ~isempty(period_type)
    sel = [sel ',period_type=' period_type];
end
measure = measures(find_record(measures,sel));
if isempty(measure)
    logmsg(['Could not find pre-cluster ' num2str(pre_cluster) ' in measures of ' recordfilter(record)])
    return
end

pairs = measure.pairs;
pair = pairs(find([pairs.post_cluster]==post_cluster,1));
if isempty(pair)
    logmsg(['Could not find post-cluster ' num2str(post_cluster) ' in measures of ' recordfilter(record)])
    return
end

ind_pre = find([sAP.sCluster.Cluster]==pre_cluster);
ind_post = find([sAP.sCluster.Cluster]==post_cluster);

tbin_centers = -0.100:params.dt:0.100;
xcortbins = [min(tbin_centers)-params.dt/2 tbin_centers+params.dt/2];
%period_type = 'opto_test';
xcorrv = zeros(length(tbin_centers),1);
duration = 0;
ind_i_all = [];
ind_j_all = [];

[spont_start,spont_end,spont_duration] = nj_get_spont_interval(record,sAP);
pulse_times = nj_get_pulse_times(record,sAP);
n_pulses = 0;
for t = record.opto_test(:)'
    if isempty(pulse_times{t})
        continue
    end
    start_time = pulse_times{t}(1) - 0.1;
    end_time = pulse_times{t}(end) + 0.1;
    duration = duration + end_time - start_time;
    n_pulses = n_pulses + length(pulse_times{t});
    ind = find(sAP.sCluster(ind_pre).SpikeTimes>start_time & sAP.sCluster(ind_pre).SpikeTimes<end_time) ;
    ind_i_all = [ind_i_all;ind];
    ind = find(sAP.sCluster(ind_post).SpikeTimes>start_time & sAP.sCluster(ind_post).SpikeTimes<end_time) ;
    ind_j_all = [ind_j_all;ind];
end
sp_i = sAP.sCluster(ind_pre).SpikeTimes(ind_i_all);
sp_j = sAP.sCluster(ind_post).SpikeTimes(ind_j_all);

if ~isempty(sp_i) &&  ~isempty(sp_j)
    xcorrv = xcorrspiketimes(sp_i,sp_j,xcortbins);
else
    xcorrv = zeros(size(tbin_centers));
end

%% Pulse-triggered rasters

n_cols = 3;
n_rows = ceil(length(pulse_times)/n_cols);
figure('Name','Pulse-triggered','NumberTitle','off')
for i=1:length(pulse_times)
    subplot(n_rows,n_cols,i)
    hold on
    rastergram(sp_i,pulse_times{i},[-0.02 0.1],Color=params.clr_l5);
    rastergram(sp_j,pulse_times{i},[-0.02 0.1],Color=params.clr_pom);
    title(sAP.cellBlock{i}.strExpType)
end


%% First spike-triggered rasters
switch period_type
    case 'best'
        sp_i = sAP.sCluster(ind_pre).SpikeTimes;
        sp_j = sAP.sCluster(ind_post).SpikeTimes;
end

[bursts_i,tonics_i,n_events_i,spikes_per_burst_i,mean_burst_isi_i] = get_bursttimes(sp_i,params);
first_sp_i = union(bursts_i,tonics_i);

[bursts_j,tonics_j,n_events_j,spikes_per_burst_j,mean_burst_isi_j] = get_bursttimes(sp_j,params);
first_sp_j = union(bursts_j,tonics_j);

if ~isempty(first_sp_i) &&  ~isempty(first_sp_j)
    first_xcorrv = xcorrspiketimes(first_sp_i,sp_j,xcortbins);
else
    first_xcorrv = zeros(size(tbin_centers));
end

% take a subset for visualization
n_subset = min(100,length(first_sp_i));
ind_subset = randi(length(first_sp_i),[n_subset 1]);
%ind_subset = 1:length(first_sp_i);

figure('Name','First spike-triggered','NumberTitle','off');
axes('units','centimeters','position',[2 1 4 1])
bar(tbin_centers,first_xcorrv)
xlim([-0.02 0.02])
ylabel('Count')
xlabel('Time from L5 spike or burst (s)')

axes('units','centimeters','position',[2 3 4 4])
hold on
[~,~,h_l5] = rastergram(sp_i,first_sp_i(ind_subset),[-0.02 0.02],"Color",params.clr_l5,'LineWidth',2);
[~,~,h_pom] = rastergram(sp_j,first_sp_i(ind_subset),[-0.02 0.02],"Color",params.clr_pom,'LineWidth',2);
set(gca,'ydir','reverse')
ylabel('L5 spike or burst # (subset)')
title([record.subject ': ' num2str(pre_cluster) '->'  num2str(post_cluster)])
legend([h_l5,h_pom],'L5','PoM','location','northwest')
% for i = 1:length(pulse_times)
%     pt = pulse_times{i}(1);
%     ind = find(sp_i>pt,1);
%     if ~isempty(ind)
%         plot(xlim,(ind-1)*[1 1])
%         text(-0.01,(ind-1),num2str(i),'verticalalignment','middle')
%     end
% end


