function results_njtestrecord( record )
%results_njtestrecord. Show results for Nora's Neuropixels data
%
%  RESULTS_NJTESTRECORD( record )
%
%  2024, Alexander Heimel
%

global measures global_record
global_record = record;
evalin('base','global measures');
evalin('base','global global_record');

params = nj_default_parameters(record);

measures = record.measures;
disp(record.subject)

show_cells(record,'thalamus');

show_cells(record,'barrel');

period_type = 'opto_test';
mask_opto = [measures.good_bc]  & contains({measures.period_type},period_type) & [measures.rate]>0;
n_pairs = sum(arrayfun(@(m) length(m.pairs), measures(mask_opto)));
disp(['Number of pairs (opto_test): ' num2str(n_pairs)] )


%% Show pair data
if params.show_pairs
    for m=1:length(measures)
        show_pairs(measures(m),params);
    end
end
%

logmsg('Measures available in workspace as ''measures'', record as ''global_record''.');
end

function show_cells(record,area)
ind = find(contains({record.measures.area},area));
measures = record.measures(ind);
figure('Name',[record.subject ' ' area],'NumberTitle','off');

mask_bad = ~[measures.good_bc] & ~[measures.good_ks];
mask_good_bc_only = [measures.good_bc] & ~[measures.good_ks];
mask_good_ks_only = [measures.good_ks] & ~[measures.good_bc];
mask_good = [measures.good_bc] & [measures.good_ks];

% Some counts

period_type = 'spont_test';
disp(['Area:' area])
mask_spont = [measures.good_bc] & ~[measures.non_somatic] & contains({measures.period_type},period_type) & [measures.rate]>0;
n_units = sum(mask_spont);
disp(['Number of units (spont_test): ' num2str(n_units)] )



measures = measures([measures.good_bc] & ~[measures.non_somatic]);
if isempty(measures)
    return
end




% Figures
subplot(2,2,1)
hold on
ind_spont_test = contains({measures.period_type},'spont_test');
ind_opto_test = contains({measures.period_type},'opto_test');
ind_opto_plus_whisker = contains({measures.period_type},'opto_plus_whisker');

% x =  [arrayfun(@(m) m.spont_test.rate,measures(ind_spont_test))];
% y = [ arrayfun(@(m) m.opto_test.rate,measures(ind_opto_test))];
x = [measures(ind_spont_test).rate];
y = [measures(ind_opto_test).rate];
scatter(x,y);
xlabel('Rate spont (sp/s)');
ylabel('Rate opto (sp/s)');
xlim([1e-4 1e2]);
ylim([1e-4 1e2]);
set(gca,'XScale','log');
set(gca,'YScale','log');
axis square
 xyline();


subplot(2,2,2)
x = [measures(ind_spont_test).burst_rate];
y = [measures(ind_opto_test).burst_rate];
 scatter(x,y);


% x =  [arrayfun(@(m) m.spont_test.burst_rate,measures)];
% y = [ arrayfun(@(m) m.opto_test.burst_rate,measures)];
%  scatter(x(mask_bad),y(mask_bad));
%  scatter(x(mask_good_bc_only),y(mask_good_bc_only),'filled','r');
%  scatter(x(mask_good_ks_only),y(mask_good_ks_only),'filled','b');
%  scatter(x(mask_good),y(mask_good),'filled', 'c');
xlabel('Burst rate spont (sp/s)');
ylabel('Burst rate opto (sp/s)');
set(gca,'XScale','log');
set(gca,'YScale','log');
xlim([1e-4 1e2]);
ylim([1e-4 1e2]);
axis square
 xyline();


subplot(2,2,3)
x = [measures(ind_spont_test).burst_rate];
y = [measures(ind_opto_plus_whisker).burst_rate];
scatter(x,y);


% x =  [arrayfun(@(m) m.opto_test.burst_rate,measures)];
% y = [ arrayfun(@(m) m.opto_plus_whisker.burst_rate,measures)];
%  scatter(x(mask_bad),y(mask_bad));
%  scatter(x(mask_good_bc_only),y(mask_good_bc_only),'filled','r');
%  scatter(x(mask_good_ks_only),y(mask_good_ks_only),'filled','b');
%  scatter(x(mask_good),y(mask_good),'filled', 'c');
xlabel('Burst rate opto (sp/s)');
ylabel('Burst rate opto plus_whisker (sp/s)');
set(gca,'XScale','log');
set(gca,'YScale','log');
xlim([1e-4 1e2]);
ylim([1e-4 1e2]);
axis square
 xyline();


 subplot(2,2,4)

 % axis off
 % set(gca,'units','points');
 % text(0,0,['n = ' num2str(length(measures))]);

 hist_frac = measures(1).isi_hist_frac;
 hist_edges = measures(1).isi_hist_edges;
 count = 0;
 for i=2:length(measures)
     if isempty(measures(i).isi_hist_frac)
         continue
     end
     hist_frac = hist_frac + measures(i).isi_hist_frac;
     count = count+1;
 end
 hist_frac = hist_frac/count;
 bar(hist_frac)
set(gca,'XTick',1:length(hist_frac))
set(gca,'Xticklabel',hist_edges(2:end))
end


function show_pairs(measures,params)

if measures.good_bc
    return
end

pairs = measures.pairs;
for p=1:length(pairs)
    pair = pairs(p);
   

    if pair.peak_snr<params.min_corr_peak_in_std
        continue
    end

    if ~contains(pair.post_area,'complex')
        continue
    end

    figure('Name',[measures.period_type ':' num2str(measures.cluster) '->'  num2str(pair.post_cluster) ] ,'NumberTitle','off');
    hold on
    h.cc_all = plot(pair.tbin_centers,pair.xcorrv,'k-','LineWidth',2);
    h.cc_first_spikes = plot(pair.tbin_centers,pair.xcorrv_first_spikes,'r--','LineWidth',2);
    h.cc_bursts = plot(pair.tbin_centers,pair.xcorrv_bursts,'b--','LineWidth',2);


    plot([0 0],ylim,'Color',[0.9 0.9 0])
    plot([-params.max_delay -params.max_delay],ylim,'Color',[0.9 0.9 0])
    plot([ params.max_delay  params.max_delay],ylim,'Color',[0.9 0.9 0])

    h.shuffled_cc = plot(pair.tbin_centers,pair.xcorrv_mean,'-','Color',[0.7 0.7 0.7]);
    plot(pair.tbin_centers,pair.xcorrv_mean - 2*pair.xcorrv_std,'--','Color',[0.7 0.7 0.7]);
    plot(pair.tbin_centers,pair.xcorrv_mean + 2*pair.xcorrv_std,'--','Color',[0.7 0.7 0.7]);
    plot(xlim,max(pair.xcorrv_mean) + params.min_corr_peak_in_std * max(pair.xcorrv_std) *[1 1],'--','Color',[0.7 0.7 0.7]);

    title(['Cluster ' num2str(measures.cluster) ' in ' measures.area ' vs ' num2str(pair.post_cluster) ' in ' pair.post_area ])
    xlabel(['Time spike of #' num2str(pair.post_cluster) ' after #' num2str(measures.cluster) ' (s) '])
    hold off

    legend([h.cc_all,h.cc_first_spikes,h.cc_bursts,h.shuffled_cc ],{'All spikes','First spike','Bursts','Shuffled'})
end % pair p

end