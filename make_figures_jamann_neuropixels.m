% make_figures_jamann_neuropixels
%
%    Script to create figures for neuropixels data of Jamann et al.
%    To re-analyse the data, load database, then experiment_db(db) and run
%    analysis in GUI. 
%
% 2024, Alexander Heimel

% To do show that
%
% - bursts in L5 are less likely to be directly followed by burst in
% POm, and compute this both for the bursts found in the entire
% population and for bursts found in putatively connected neuron?
%
% - if the number of spikes per burst is smaller in POm, for 1. all
% bursts, 2. bursts following L5 bursts, 3, bursts in neurons putatively
% receiving L5 input. Perhaps actually Nora already looked at 1 when she
% made figure 4g and 4h.
%
% - run Nora's original script to find the L5 cells that I am missing. In
% Nora's Fig 4d she has 192 neurons from 12 mice. I can only find 85
% (somatic) or 104 (including non-somatic) from 12 mice.
%
%
% - check for differences in ISI histograms for Ctrl and Cpz
%       -> cannot find any differences that are consistent over selectia
%          criteria, but made some reviewer figures
%
% - find pairs separately in the different trial types and then use all pairs
%       -> variability increased, but slightly increased number of pairs
% 
% - compute post-synaptic spike probability in first and second window, 
%   and see if that shifts towards the second window for the Cpz model
%       ->  use 10 ms and 20 ms windows, 
% 
% - correlate connectivity-likelihood with burst transfer
%       -> no correlation
%
% - check the low burst transfer pairs. Are these also low likelihood of connectivity
% 
% - correlate presynaptic burst frequency with burst transfer success
%        -> quite a strong correlation, but not sure what to make of it
% 
%  From discussion with Maarten
% - check if post-synaptic spike time variability is larger in Cpz
%     -> no, Nested t-test P = 0.63. n = 35 neurons/pairs, N = 5 mice (Ctrl), n = 5 neurons/pairs, N = 3 mice (Cpz).
%             (post_spike_time_std)
%
% - move burst fraction figure (with significant result for L5) to main figure
%     -> done, 'bursting neurons' here defined as all neurons with at least
%        one burst
%
% - correlate chance of spike/burst transfer for with the frequency of
%        the presynapic burst (and see if is different in Cpz animals)
%            -> burst frequency for successfully transferred is
%            significantly higher (stats per cell), than those of
%            not-transferred burst. However, relationship looks similar
%            for Ctrl and Cpz
%
% From Nora:
% -ISI distribution:
%   I don't really understand the reason for using a cutoff for the ISI fraction, 
%   rather than looking at the entire distribution and then doing e.g. a Kolmogorov Smirnow?
%   Also did you mean 1ms (100 Hz, bursting) in the y axis? Or why did you look for ISIs below 1s?
%   I would love to see how the distribution looks and what the differences are.
%
% - Could you share the figures as vector graphics? Or even better share the 
%   raw data in the graphs, so I can make them in the same GraphPad style 
%   that I used for the other figures (although its already really close!)?
%
% Still to do (following discussion with Maarten 2025-03-31)
%
% - email additional figure about changed burst frequency in L5
%
% - in the paired data:
%   split bursts with frequency above 200 Hz, and those below.
%   compute transfer rates for these two groups separately
%   see if fast bursts are more successful, and perhaps also the fraction
%   of fast bursts
%   Hypothesis: in Cpz model, the transfer rates of the fast bursts have
%   gone down to the level of the slow burst in Ctrl

CTRL = 1;
CPZ = 2;
conditions = {'Ctrl','Cpz'};
verbose = true;
savefigs = true;

params = nj_default_parameters();
filename = fullfile(params.projectfolder,'Data_collection','NPXdataAlexander','nj_db.mat');
load(filename,'db');
disp(['Loaded ' filename])

if params.stringent_cell_quality
    add_crit = ',good_bc=1,non_somatic=0,rate>0';
else
    add_crit = '';
end
logmsg(['params.stringent_cell_quality = ' num2str(params.stringent_cell_quality)]);

return

%% Open experiment db
experiment_db(filename);

%% Duplicate Figure 4d - L5 Firing frequency Ctrl vs Cpz 
tbl = nj_collect_values( db, conditions,{'rate'},['period_type=spont_test,area=*barrel*' add_crit]);
nj_all_data_comparison(tbl,'rate','Firing frequency (Hz)','log',[1e-4 inf]); %#ok<*UNRCH>
title('L5','FontWeight','normal','Color',params.clr_l5);
if savefigs
    exportgraphics(gcf,'fig4d.pdf');
end

%% Duplicate Figure 4e - POm Firing frequency Ctrl vs Cpz 
tbl = nj_collect_values( db, conditions,'rate',['period_type=spont_test,area=*complex*' add_crit]);
nj_all_data_comparison(tbl,'rate','Firing frequency (Hz)','log',[1e-4 inf]); %#ok<*UNRCH>
title('POm','FontWeight','normal','Color',params.clr_pom)
if savefigs
    exportgraphics(gcf,'fig4e.pdf');
end

%% Fig 4g - L5 Fraction of burst Ctrl vs Cpz 
tbl = nj_collect_values( db, conditions,{'burst_fraction'},['period_type=spont_test,area=*barrel*' add_crit]);
tbl(tbl.burst_fraction==0,:) = []; % only consider neurons that have at least 1 burst
nj_all_data_comparison(tbl,'burst_fraction','Fraction of bursts'); %#ok<*UNRCH>
title('L5','FontWeight','normal','Color',params.clr_l5);
if savefigs
     exportgraphics(gcf,fullfile(getdesktopfolder(),'fig4g.pdf'));
     writetable(tbl,fullfile(getdesktopfolder(),'fig4g.xlsx'))
end

%% Fig 4h - POm Fraction of burst Ctrl vs Cpz 
tbl = nj_collect_values( db, conditions,{'burst_fraction','bursting'},['period_type=spont_test,area=*complex*' add_crit]);
%tbl(~tbl.bursting,:) = [];
tbl(tbl.burst_fraction==0,:) = []; % only consider neurons that have at least 1 burst
nj_all_data_comparison(tbl,'burst_fraction','Fraction of bursts'); %#ok<*UNRCH>
title('POm','FontWeight','normal','Color',params.clr_pom);
if savefigs
     exportgraphics(gcf,fullfile(getdesktopfolder,'fig4h.pdf'));
     writetable(tbl,fullfile(getdesktopfolder(),'fig4h.xlsx'))
end

%% Duplicate Figure 4i - POm Bursting fraction Ctrl vs Cpz 
tbl = nj_collect_values( db, conditions,'bursting',['period_type=spont_test,area=*complex*' add_crit]);
nj_all_data_comparison(tbl,'bursting','Fraction of bursting units',[],[],[],0); %#ok<*UNRCH>
title('POm','FontWeight','normal','Color',params.clr_pom);
if savefigs
    exportgraphics(gcf,'fig4g.pdf');
end

%% Duplicate Figure 4j - POm Burst frequency Ctrl vs Cpz 
tbl = nj_collect_values( db, conditions,{'burst_frequency','bursting'},['period_type=spont_test,area=*complex*,bursting=1' add_crit]);
tbl(~tbl.bursting,:) = [];
nj_all_data_comparison(tbl,'burst_frequency','Burst frequency (Hz)','linear',[0 500]); %#ok<*UNRCH>
title('POm','FontWeight','normal','Color',params.clr_pom);
if savefigs
    exportgraphics(gcf,'fig4j.pdf');
end

%% Extra figure - L5 Burst frequency Ctrl vs Cpz 
tbl = nj_collect_values( db, conditions,{'burst_frequency','bursting'},['period_type=spont_test,area=*barrel*' add_crit]);
%tbl(~tbl.bursting,:) = [];
nj_all_data_comparison(tbl,'burst_frequency','Burst frequency (Hz)','linear',[0 500]); %#ok<*UNRCH>
title('L5','FontWeight','normal','Color',params.clr_l5);
if savefigs
    exportgraphics(gcf,'fig_burst_frequency_l5.pdf');
    writetable(tbl,fullfile(getdesktopfolder(),'fig_burst_frequency_l5.xlsx'))
end

%% Extra Figure 4x - POm Burst length Ctrl vs Cpz 
tbl = nj_collect_values( db, conditions,'n_spikes_per_burst',['period_type=opto_test,area=*complex*,bursting=1' add_crit]);
nj_all_data_comparison(tbl,'n_spikes_per_burst','Spikes per burst','linear',[0 500]); %#ok<*UNRCH>
title('POm','FontWeight','normal','Color',params.clr_pom)
if savefigs
    exportgraphics(gcf,'fig_spikes_per_burst_pom.pdf');
    writetable(tbl,fullfile(getdesktopfolder(),'fig_spikes_per_burst_pom.xlsx'))
end

%% Find examples for Figure 4k
%tbl = nj_collect_pair_values(db,conditions,{'peak_time','peak_snr','post_cluster','post_cluster_good_bc','post_area'},{'area','cluster'},'area=*layer*,period_type=opto_test,good_bc=1');
tbl = nj_collect_pair_values(db,conditions,{'peak_time','peak_snr','peak_height','post_cluster','post_cluster_good_bc','post_area'},{'area','cluster'},'area=*layer*,period_type=best,good_bc=1');


tbl(~tbl.post_cluster_good_bc,:) = [];
tbl(~contains(tbl.post_area,'complex'),:) = [];
tbl(tbl.peak_snr<params.min_corr_peak_in_std,:) = [];
tbl.cluster = cell2mat(tbl.cluster);

% find Ctrl examples
tbl(~(tbl.condition=="Ctrl"),:) = [];
%tbl = tbl(tbl.peak_time<0.011 & tbl.peak_time>0.004,:);
%tbl = tbl(tbl.peak_snr>6,:);
tbl = tbl(tbl.peak_height>20 & tbl.peak_height<=40,:);

% find Cpz examples
% tbl(~(tbl.condition=="Cpz"),:) = [];
% tbl = tbl(tbl.peak_time>0.011,:);
% tbl = tbl(tbl.peak_snr>6,:);

%% Example finder
c = 4
nj_plot_pair_rasters(db(find_record(db,['subject=' char(tbl.subject(c))])),tbl.cluster(c),tbl.post_cluster(c),'best') % 


%% Extra figure 4k Example connected pairs rastergram
% Ctrl: 1,4,12 - 4: Iv1x4 1139 -> 711
record = db(find_record(db,'subject=Iv1x4'));
measures = record.measures;
measure = measures(find_record(measures,'cluster=1139,period_type=opto_test'));
pairs = measure.pairs;
pair = pairs(find([pairs.post_cluster]==711));
sAP = nj_load_data(record);

%% Ctrl example raster
% nj_plot_pair_rasters(db(find_record(db,'subject=Iv3a8')),1277,76,'opto_test') % test 1 looks good (10 spikes)
% nj_plot_pair_rasters(db(find_record(db,'subject=Iv4a11')),272,13,'opto_test') % peak 9 spikes, some other clusters post to 272 are probably the same unit
% nj_plot_pair_rasters(db(find_record(db,'subject=Iv1x4')),1139,711,'opto_test') % first example
% nj_plot_pair_rasters(db(find_record(db,'subject=Iv3a8')),1598,76,'opto_test')
% nj_plot_pair_rasters(db(find_record(db,'subject=Iv3a8')),1687,1012,'best') % not opto, many spikes
% nj_plot_pair_rasters(db(find_record(db,'subject=Iv1x10')),357,306,'best') % not opto, many spikes
% nj_plot_pair_rasters(db(find_record(db,'subject=Iv3a8')),1277,76,'best') % not opto
% % unit 76 and 1034 are probably coming from the same unit
% nj_plot_pair_rasters(db(find_record(db,'subject=Iv3a8')),1598,1034,'opto_test') % 
nj_plot_pair_rasters(db(find_record(db,'subject=Iv3a8')),1687,1012,'best')

%% Cpz example raster
% nj_plot_pair_rasters(db(find_record(db,'subject=Iv2b7')),582,16,'opto_test'); % alternative example
nj_plot_pair_rasters(db(find_record(db,'subject=Iv2b7')),582,442,'opto_test'); % looks nice. jitter between pre and post is smaller than between pulse and pre

%% Figure 4k Example connected pairs
% Ctrl
record = db(find_record(db,'subject=Iv3a8'));
measures = record.measures;
measure = measures(find_record(measures,'cluster=1687,period_type=best'));
pairs = measure.pairs;
pair = pairs(find([pairs.post_cluster]==1012));
figure;
subplot('position',[0 0.5 0.5 0.5]);
hold on
bar(pair.tbin_centers,pair.xcorrv_mean + pair.xcorrv_std,...
    'FaceColor',[0.7 0.7 0.7],'BarWidth',1,'LineStyle','none','basevalue',0);
bar(pair.tbin_centers,pair.xcorrv_mean - pair.xcorrv_std,...
    'FaceColor',[1.0 1.0 1.0],'BarWidth',1,'LineStyle','none','BaseValue',0);
stairs(pair.tbin_centers-(pair.tbin_centers(2)-pair.tbin_centers(1))/2,pair.xcorrv_mean,'Color',[0.4 0.4 0.4]);
b = bar(pair.tbin_centers,pair.xcorrv_first_spikes,'k','BarWidth',1,'BaseValue',0)
yl = ylim();
ylim([0 yl(2)])
xlim([-0.02 0.02])
set(gca,'xtick',[])
ylabel('Count');
xl = xlim;
text(xl(1),yl(2),'  Ctrl','Color',params.clr_ctrl,'FontWeight','normal','VerticalAlignment','top',HorizontalAlignment='left');
set(gca,'units','centimeters','position',[3 4. 4 1.5]);

% Cpz
record = db(find_record(db,'subject=Iv2b7'));
measures = record.measures;
measure = measures(find_record(measures,'cluster=582,period_type=opto_test'));
pairs = measure.pairs;
pair = pairs(find([pairs.post_cluster]==442));
subplot('position',[0 0 0.5 0.2]);

hold on
bar(pair.tbin_centers,pair.xcorrv_mean + pair.xcorrv_std,...
    'FaceColor',[0.7 0.7 0.7],'BarWidth',1,'LineStyle','none','basevalue',0);
bar(pair.tbin_centers,pair.xcorrv_mean - pair.xcorrv_std,...
    'FaceColor',[1.0 1.0 1.0],'BarWidth',1,'LineStyle','none','BaseValue',0);
stairs(pair.tbin_centers-(pair.tbin_centers(2)-pair.tbin_centers(1))/2,pair.xcorrv_mean,'Color',[0.4 0.4 0.4]);
b = bar(pair.tbin_centers,pair.xcorrv_first_spikes,'BarWidth',1,'BaseValue',0,'FaceColor',params.clr_cpz, 'EdgeColor',params.clr_cpz)
yl = ylim();
ylim([0 yl(2)])
xlim([-0.02 0.02])
xlabel('Time from L5 to POm spike (s)');
ylabel('Count');
yl = ylim;
xl = xlim;
text(xl(1),yl(2),'  Cpz','Color',params.clr_cpz,'FontWeight','normal','VerticalAlignment','top',HorizontalAlignment='left');
set(gca,'units','centimeters','position',[3 2 4 1.5]);

if savefigs
    exportgraphics(gcf,fullfile(getdesktopfolder(),'fig4k.pdf'));
end

%% Figure 4l Connection probability
connection_probability = {};
for c = 1:2
    ind = find(arrayfun(@(x) strcmp(x.condition,conditions{c}),db));
    n_mice = length(ind);
    n_l5_cells{c} = zeros(n_mice,1);
    n_pom_cells{c} = zeros(n_mice,1);
    n_potential_connections{c} = zeros(n_mice,1);
    n_connections{c} = zeros(n_mice,1);
    connection_probability{c} = NaN(n_mice,1);
    
    for i = 1:length(ind)
        measures = db(ind(i)).measures;

        if params.stringent_cell_quality
            measures = measures(find_record(measures,'good_bc=1'));
        end

        if params.stringent_post_area
            n_pom_cells{c}(i) = sum(contains({measures.area},'complex'));
        else
            n_pom_cells{c}(i) = sum(contains({measures.area},'thalamus'));
        end

        n_l5_cells{c}(i) = sum(contains({measures.area},'layer 5'));
        n_potential_connections{c}(i) =  n_l5_cells{c}(i) * n_pom_cells{c}(i);
        ind_l5 = find(contains({measures.area},'layer 5'));
        ind_l5 = find_record(measures,'area=*layer 5*,period_type=best');
        n_connections{c}(i) = 0;
        for j = 1:length(ind_l5)
            pairs = measures(ind_l5(j)).pairs;
            if params.stringent_post_area
                pairs = pairs(find_record(pairs,'post_area=*complex*'));
            end
            if params.stringent_cell_quality
                pairs = pairs(find_record(pairs,'post_cluster_good_bc=1'));
            end
            pairs = pairs(find_record(pairs,['peak_snr>' num2str(params.min_corr_peak_in_std)]));

            n_connections{c}(i) = n_connections{c}(i) + length(pairs);
        end
        if n_potential_connections{c}(i)>0
            connection_probability{c}(i) = n_connections{c}(i) /  n_potential_connections{c}(i) ;
        end
    end
end
connection_probability_perc = cellfun(@(x) x*100,connection_probability,'UniformOutput',false);
h = ivt_graph(connection_probability_perc,[],'ylab','Connection prob. (%)',...
    'xticklabels',{'Ctrl','Cpz'},'spaced',3,'color',{[1 1 1],[0 0 0]},...
    'style','level','errorbars_sides','both','markersize',6,'barwidth',0.5,...
    'test','ttest2');
set(gca,'units','centimeters','position',[2 1 1.5 3.5])
ylim([0 0.2])
fontsize(scale=0.8);
n_ctrl = sum(~isnan(connection_probability{CTRL}));
n_cpz = sum(~isnan(connection_probability{CPZ}));


disp(['Unpaired t-test P = ' num2str(h.p_sig{2},2) ...
    '. N = ' num2str(n_ctrl) ' mice (Ctrl), N = ' num2str(n_cpz) ' mice (Cpz).'] )

condition = [repmat("Ctrl",n_ctrl,1);repmat("Cpz",n_cpz,1)]
connection_probability = [...
    connection_probability_perc{1}(~isnan(connection_probability_perc{1})) ; ...
    connection_probability_perc{2}(~isnan(connection_probability_perc{2}))]
tbl = table(condition,connection_probability)

if savefigs
    exportgraphics(gcf,fullfile(getdesktopfolder(),'figl.pdf'));
    writetable(tbl,fullfile(getdesktopfolder(),'fig4l.xlsx'))
end

%% Alternative Figure 4o Delay - to POM
params = nj_default_parameters();
tbl = nj_collect_pair_values(db,conditions,{'peak_time','peak_snr','peak_height','post_cluster_good_bc','post_area'},[],'area=*layer*,period_type=best,good_bc=1');
tbl(~contains(tbl.post_area,'complex') ,:) = [];
tbl(tbl.peak_snr<params.min_corr_peak_in_std,:) = [];
tbl.peak_time = tbl.peak_time * 1000; % ms
nj_all_data_comparison(tbl,'peak_time','Delay (ms)',[],[],[0 20]); %#ok<*UNRCH>
title('L5\rightarrowPOm','FontWeight','normal');
if savefigs
    exportgraphics(gcf,fullfile(getdesktopfolder(),'fig4m.pdf'));
    writetable(tbl,fullfile(getdesktopfolder(),'fig4m.xlsx'))
end

%% Old Figure 4n Fraction of spikes of first 20 ms post pre-synaptic spikes, that occur in first 10 ms.
params = nj_default_parameters();
tbl = nj_collect_pair_values( db, conditions,{'fraction_post_10ms','peak_snr','post_cluster_good_bc','post_area'},[],'area=*layer*,period_type=best,good_bc=1');
if params.stringent_post_area 
    tbl(~contains(tbl.post_area,'complex'),:) = [];
end
tbl(tbl.peak_snr<params.min_corr_peak_in_std,:) = [];
nj_all_data_comparison(tbl,'fraction_post_10ms','Fraction in first 10 ms'); %#ok<*UNRCH>
if savefigs
    exportgraphics(gcf,fullfile(getdesktopfolder(),'fig4n.pdf'));
    writetable(tbl,fullfile(getdesktopfolder(),'fig4n.xlsx'))
end

%% Alternative Figure 4n Fraction of spikes of first 20 ms post pre-synaptic spikes, that occur in first 10 ms after the first spike 
params = nj_default_parameters();
tbl = nj_collect_pair_values( db, conditions,{'fraction_10ms_post_first','peak_snr','post_cluster_good_bc','post_area'},[],'area=*layer*,period_type=best,good_bc=1');
if params.stringent_post_area 
    tbl(~contains(tbl.post_area,'complex'),:) = [];
end
tbl(tbl.peak_snr<params.min_corr_peak_in_std,:) = [];
nj_all_data_comparison(tbl,'fraction_10ms_post_first','Fraction in first 10 ms after first spike'); %#ok<*UNRCH>
if savefigs
    exportgraphics(gcf,fullfile(getdesktopfolder(),'fig4n_10ms_post_first_spike.pdf'));
    writetable(tbl,fullfile(getdesktopfolder(),'fig4n_10ms_post_first_spike.xlsx'))
end

%% Figure 4o_delay_first
params = nj_default_parameters();
tbl = nj_collect_pair_values( db, conditions,{'delay_first_increased_spikes','peak_snr','post_cluster_good_bc','post_area'},[],'area=*layer*,period_type=best,good_bc=1');
if params.stringent_post_area 
    tbl(~contains(tbl.post_area,'complex'),:) = [];
end
tbl(tbl.peak_snr<params.min_corr_peak_in_std,:) = [];
tbl.delay_first_increased_spikes = tbl.delay_first_increased_spikes * 1000;
nj_all_data_comparison(tbl,'delay_first_increased_spikes','Delay first spikes (ms)'); %#ok<*UNRCH>
ylim([0 15]);
if savefigs
    exportgraphics(gcf,fullfile(getdesktopfolder(),'fig4x_delay_first.pdf'));
    writetable(tbl,fullfile(getdesktopfolder(),'fig4x_delay_first.xlsx'))
end

%% Extra Figure 4x std in delay first spike
params = nj_default_parameters();
tbl = nj_collect_pair_values( db, conditions,{'post_first_spike_time_std','peak_snr','post_cluster_good_bc','post_area'},[],'area=*layer*,period_type=best,good_bc=1');
if params.stringent_post_area 
    tbl(~contains(tbl.post_area,'complex'),:) = [];
end
tbl(tbl.peak_snr<params.min_corr_peak_in_std,:) = [];
tbl.post_first_spike_time_std = tbl.post_first_spike_time_std * 1000;
nj_all_data_comparison(tbl,'post_first_spike_time_std','Std. dev. in first spike (ms)'); %#ok<*UNRCH>
%ylim([0 15]);
if savefigs
    exportgraphics(gcf,fullfile(getdesktopfolder(),'fig4x_std_first.pdf'));
    writetable(tbl,fullfile(getdesktopfolder(),'fig4x_std_first.xlsx'))
end

%% Extra Figure 4x variance in delay first spike
params = nj_default_parameters();
tbl = nj_collect_pair_values( db, conditions,{'post_first_spike_time_var','peak_snr','post_cluster_good_bc','post_area'},[],'area=*layer*,period_type=best,good_bc=1');
if params.stringent_post_area 
    tbl(~contains(tbl.post_area,'complex'),:) = [];
end
tbl(tbl.peak_snr<params.min_corr_peak_in_std,:) = [];
tbl.post_first_spike_time_var = tbl.post_first_spike_time_var * 1000^2;
nj_all_data_comparison(tbl,'post_first_spike_time_var','Variance in first spike (ms^2)'); %#ok<*UNRCH>
%ylim([0 15]);
if savefigs
    exportgraphics(gcf,fullfile(getdesktopfolder(),'fig4x_var_first.pdf'));
    writetable(tbl,fullfile(getdesktopfolder(),'fig4x_var_first.xlsx'))
end


%% Extra Figure 4x Burst transfer
params = nj_default_parameters();
tbl = nj_collect_pair_values( db, conditions,{'burst_transfer','peak_snr','post_cluster_good_bc','post_area'},[],'area=*layer*,period_type=best,good_bc=1');
if params.stringent_post_area 
    tbl(~contains(tbl.post_area,'complex'),:) = [];
end
tbl(tbl.peak_snr<params.min_corr_peak_in_std,:) = [];

% logit transform
tbl.burst_transfer = tbl.burst_transfer+0.01; % there is a lower bound on the burst_transfer that we can measure, ideally should be taken from data
tbl.burst_transfer = log(tbl.burst_transfer./(1-tbl.burst_transfer));

nj_all_data_comparison(tbl,'burst_transfer','Burst transfer'); %#ok<*UNRCH>
if savefigs
    exportgraphics(gcf,fullfile(getdesktopfolder(),'fig4x_burst_transfer.pdf'));
    writetable(tbl,fullfile(getdesktopfolder(),'fig4x_burst_transfer.xlsx'))
end


%% Extra Figure 4x Spikes post burst
params = nj_default_parameters();
tbl = nj_collect_pair_values( db, conditions,{'spikes_post_burst','peak_snr','post_cluster_good_bc','post_area'},[],'area=*layer*,period_type=best,good_bc=1');
if params.stringent_post_area 
    tbl(~contains(tbl.post_area,'complex'),:) = [];
end
tbl(tbl.peak_snr<params.min_corr_peak_in_std,:) = [];
nj_all_data_comparison(tbl,'spikes_post_burst','spikes_post_burst'); %#ok<*UNRCH>
if savefigs
    exportgraphics(gcf,fullfile(getdesktopfolder(),'fig4x_spikes_post_burst.pdf'));
    writetable(tbl,fullfile(getdesktopfolder(),'fig4x_spikes_post_burst.xlsx'))
end

%% Extra Figure 4x Spike transfer
params = nj_default_parameters();
tbl = nj_collect_pair_values( db, conditions,{'spike_transfer','peak_snr','post_cluster_good_bc','post_area'},[],'area=*layer*,period_type=best,good_bc=1');
if params.stringent_post_area 
    tbl(~contains(tbl.post_area,'complex'),:) = [];
end
tbl(tbl.peak_snr<params.min_corr_peak_in_std,:) = [];
nj_all_data_comparison(tbl,'spike_transfer','spike_transfer'); %#ok<*UNRCH>
if savefigs
    exportgraphics(gcf,fullfile(getdesktopfolder(),'fig4x_spike_transfer.pdf'));
    writetable(tbl,fullfile(getdesktopfolder(),'fig4x_spike_transfer.xlsx'))
end


%% Figure 4x Burst freq of transferred vs not-transferred spikes
params = nj_default_parameters();
tbl = nj_collect_pair_values( db, conditions,{'transferred_burst_isi','not_transferred_burst_isi','peak_snr','post_cluster_good_bc','post_area'},[],'area=*layer*,period_type=best,good_bc=1');
if params.stringent_post_area 
    tbl(~contains(tbl.post_area,'complex'),:) = [];
end
tbl(tbl.peak_snr<params.min_corr_peak_in_std,:) = [];

figure;
hold on
scatter(1./tbl.not_transferred_burst_isi(tbl.condition=="Ctrl"),...
    1./tbl.transferred_burst_isi(tbl.condition=="Ctrl"),...
    'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[1 1 1]*0.6)
scatter(1./tbl.not_transferred_burst_isi(tbl.condition=="Cpz"),...
    1./tbl.transferred_burst_isi(tbl.condition=="Cpz"),...
    'MarkerFaceColor',[1 1 1]*0.6,'MarkerEdgeColor',[1 1 1]*0.6)
axis square
xlabel('Freq. of not transferred burst')
ylabel('Freq. of transferred burst')
ylim([100 400]);
xlim([100 400]);
xyline
[h,p] = ttest(1./tbl.not_transferred_burst_isi,1./tbl.transferred_burst_isi)
plot_significance(250,250,380,p)
if savefigs
     exportgraphics(gcf,fullfile(getdesktopfolder,'fig4x_burst_freq_of_success.pdf'));
end

%% Figure 4x Transfer fast bursts relative to slow burst
params = nj_default_parameters();
tbl = nj_collect_pair_values( db, conditions,{'fast_burst_transfer','slow_burst_transfer','peak_snr','post_cluster_good_bc','post_area'},[],'area=*layer*,period_type=best,good_bc=1');
if params.stringent_post_area 
    tbl(~contains(tbl.post_area,'complex'),:) = [];
end
tbl(tbl.peak_snr<params.min_corr_peak_in_std,:) = [];
tbl.relative_transfer = tbl.fast_burst_transfer - tbl.slow_burst_transfer;
nj_all_data_comparison(tbl,'relative_transfer','Relative transfer'); %#ok<*UNRCH>

%% Non-hierarchical non-parametric stats
x = tbl.relative_transfer(tbl.condition=="Ctrl");
y = tbl.relative_transfer(tbl.condition=="Cpz");
p = ranksum(x,y);
disp(['Difference between Ctrl and Cpz, rank-sum p = ' num2str(p,2)])
p = signrank(x);
disp(['Effect fast vs slow for Ctrl, signrank p = ' num2str(p,2)])
pwr = bootstrap_power(x,length(y),'ttest',0.05);  % actually, posthoc_poweranalysis is more generic
disp(['Power for same test with n Cpz pairs = ' num2str(pwr,2)])
pwr = bootstrap_power(x,32,'ttest',0.05)

%% Figure 4x Burst transferred
params = nj_default_parameters();
tbl = nj_collect_pair_values( db, conditions,{'fast_burst_transfer','slow_burst_transfer','peak_snr','post_cluster_good_bc','post_area'},[],'area=*layer*,period_type=best,good_bc=1');
if params.stringent_post_area 
    tbl(~contains(tbl.post_area,'complex'),:) = [];
end
tbl(tbl.peak_snr<params.min_corr_peak_in_std,:) = [];
tbl.fast_burst_transfer = 100 * tbl.fast_burst_transfer; 
tbl.slow_burst_transfer = 100 * tbl.slow_burst_transfer; 

figure
tbl(isnan(tbl.fast_burst_transfer),:) = [];
tbl_per_subject = groupsummary(tbl,{'subject', 'condition'},"mean",{'fast_burst_transfer','slow_burst_transfer'});
tbl_per_condition = groupsummary(tbl_per_subject,'condition',"mean",{'mean_fast_burst_transfer','mean_slow_burst_transfer'});
ax = axes();
ivt_graph(...
    {tbl.fast_burst_transfer(tbl.condition=="Ctrl"),...
    tbl.slow_burst_transfer(tbl.condition=="Ctrl")},[1 2],...
    'style','level','spaced',3,'test','none','markersize',3,'barwidth',0,...
    'errorbars','none','markers','open_circle','color',[1 0.7 0.7],...
    'showpairing',true,...
    'axishandle',ax);

ivt_graph(...
    {tbl.fast_burst_transfer(tbl.condition=="Cpz"),...
    tbl.slow_burst_transfer(tbl.condition=="Cpz")},[4 5],...
    'style','level','spaced',3,'test','none','markersize',3,'barwidth',0,...
    'errorbars','none','markers','open_circle','color',[0 1 0],...
    'showpairing',true,...
    'axishandle',ax);

c=get(gca,'children');
for i=1:length(c)
    set(c(i),'linewidth',0.5)
    set(c(i),'color',0.7*[1 1 1])
end

% means
ivt_graph({...
    tbl_per_subject.mean_fast_burst_transfer(tbl_per_subject.condition=="Ctrl"),...
    tbl_per_subject.mean_slow_burst_transfer(tbl_per_subject.condition=="Ctrl")},[1 2],...
    'color',{[1 1 1],[1 1 1]},'test','none','errorbars','none',...
    'style','level','spaced',3,'axishandle',gca,'markersize',6,'barwidth',0);

ivt_graph({...
    tbl_per_subject.mean_fast_burst_transfer(tbl_per_subject.condition=="Cpz"),...
    tbl_per_subject.mean_slow_burst_transfer(tbl_per_subject.condition=="Cpz")},[4 5],...
    'color',{[0 0 0],[0 0 0]},'test','none','errorbars','none',...
    'style','level','spaced',3,'axishandle',gca,'markersize',6,'barwidth',0);

% mean of means

% all conditions
ivt_graph({...
    tbl_per_condition.mean_mean_fast_burst_transfer(tbl_per_condition.condition=="Ctrl"),...
    tbl_per_condition.mean_mean_slow_burst_transfer(tbl_per_condition.condition=="Ctrl")},[1 2],...
    'color',{[1 1 1],[1 1 1]},'test','none',...
    'style','level','spaced',3,'axishandle',gca,'markersize',6,...
    'barwidth',0.5,'showpoints',0);

ivt_graph({...
    tbl_per_condition.mean_mean_fast_burst_transfer(tbl_per_condition.condition=="Cpz"),...
    tbl_per_condition.mean_mean_slow_burst_transfer(tbl_per_condition.condition=="Cpz")},[4 5],...
    'color',{[0 0 0],[0 0 0]},'test','none',...
    'style','level','spaced',3,'axishandle',gca,'markersize',6,...
    'barwidth',0.5,'showpoints',0);

ylabel('Bursts transferred (%)')
ylim([0 140])
xlim([0 6])
set(gca,'ytick',(0:20:100))
set(gca,'units','centimeters','position',[3 2 3.5 4.0]);
set(gca,'xtick',[1 2 4 5]);
set(gca,'xticklabel',{'Fast','Slow','Fast','Slow'})
set(gca,'XTickLabelRotation',60)
text(1.5,140,'Ctrl','horizontalalignment','center');
text(4.5,140,'Cpz','horizontalalignment','center');

% Statistics
tbl2 = table();
tbl2.subject = [tbl.subject;tbl.subject];
tbl2.pair_id = [ (1:height(tbl))'; (1:height(tbl))'];
tbl2.condition = [tbl.condition;tbl.condition];
tbl2.transfer = [tbl.fast_burst_transfer;tbl.slow_burst_transfer];
tbl2.burst_speed = [repmat("fast",height(tbl),1);repmat("slow",height(tbl),1)];

lme = fitlme(tbl2(tbl2.condition=="Ctrl",:),'transfer~burst_speed+ (1|subject) + (1|pair_id)');
[~,~,fe] = fixedEffects(lme);
p = fe.pValue(string(fe.Name)=="burst_speed_slow");
disp(['Bursts transferred fast vs slow in Ctrl, Hierarchical t-test: p = ' num2str(p,2) ', '])
plot_significance(1,2,110,p)

lme = fitlme(tbl2(tbl2.condition=="Cpz",:),'transfer~burst_speed+ (1|subject) + (1|pair_id)');
[~,~,fe] = fixedEffects(lme);
p = fe.pValue(string(fe.Name)=="burst_speed_slow");
disp(['Bursts transferred fast vs slow in Cpz, Hierarchical t-test: p = ' num2str(p,2) ', '])
plot_significance(4,5,110,p,[],[],[],'ns')

lme = fitlme(tbl2,'transfer~condition+burst_speed+burst_speed:condition + (1|subject) + (1|pair_id)');
[~,~,fe] = fixedEffects(lme);
fprintf('Full LME: ')
fprintf(['Ctrl vz Cpz: p = ' num2str(fe.pValue(string(fe.Name)=="condition_Ctrl"),2) ', '])
fprintf(['Fast vs Slow: p = ' num2str(fe.pValue(string(fe.Name)=="burst_speed_slow"),'%.2f') ', '])
fprintf(['Interaction: p = ' num2str(fe.pValue(string(fe.Name)=="condition_Ctrl:burst_speed_slow"),2)])
fprintf('\n')

if savefigs
    exportgraphics(gcf,fullfile(getdesktopfolder(),'fig_burst_transferred.pdf'));
    writetable(tbl2,fullfile(getdesktopfolder(),'fig_burst_transferred.xlsx'))
end


%% Looking at difference
tbl.diff = tbl.fast_burst_transfer-tbl.slow_burst_transfer;
tbl.pair_id = (1:height(tbl))';
lme = fitlme(tbl,'diff~condition + (1|subject) + (1|pair_id)');

%% Figure 4x Std in postsynaptic spike times
params = nj_default_parameters();
tbl = nj_collect_pair_values( db, conditions,{'post_spike_time_std','peak_snr','post_cluster_good_bc','post_area'},[],'area=*layer*,period_type=best,good_bc=1');
if params.stringent_post_area 
    tbl(~contains(tbl.post_area,'complex'),:) = [];
end
tbl(tbl.peak_snr<params.min_corr_peak_in_std,:) = [];
tbl.post_spike_time_std = tbl.post_spike_time_std* 1000;
nj_all_data_comparison(tbl,'post_spike_time_std','Std. dev. in delay (ms)'); %#ok<*UNRCH>
if savefigs
     exportgraphics(gcf,fullfile(getdesktopfolder,'fig4x_post_spike_time_std.pdf'));
     writetable(tbl,fullfile(getdesktopfolder(),'fig4x_post_spike_time_std.xlsx'))

end

%% Reviewer Figure 1a ISI Histogram, L5
[hist_ctrl,hist_edges,frac_below_1s_ctrl] = collect_isi(db,'Ctrl','spont_test','barrel',params);
[hist_cpz,hist_edges,frac_below_1s_cpz] = collect_isi(db,'Cpz','spont_test','barrel',params);
 figure
 hold on
 bar(hist_ctrl,'FaceColor',[1 1 1])
 bar(hist_cpz,'FaceColor',[0 0 0],'FaceAlpha',0.5)
set(gca,'XTick',(1:3:length(hist_ctrl))+0.5)
set(gca,'Xticklabel',hist_edges(2:3:end))
xlim([1.4 length(hist_ctrl)+1])
xlabel('ISI limit (s)')
ylabel('Fraction of ISIs')
title('L5','FontWeight','normal','Color',params.clr_l5)
legend('Ctrl','Cpz','Location','NorthWest')
legend boxoff
set(gca,'units','centimeters','position',[3 2 4 3.5]);
if savefigs
     exportgraphics(gcf,fullfile(getdesktopfolder,'revfig.pdf'));
end


%% Reviewer Figure 1b - L5 Fraction ISI below 1s Ctrl vs Cpz 
tbl = nj_collect_values( db, conditions,{'frac_below_1s','bursting','rate'},['period_type=spont_test,area=*barrel*' add_crit]);
tbl(tbl.frac_below_1s==0,:) = [];
nj_all_data_comparison(tbl,'frac_below_1s','Fraction of ISI below 1s'); %#ok<*UNRCH>
title('L5','FontWeight','normal','Color',params.clr_l5);
if savefigs
     exportgraphics(gcf,fullfile(getdesktopfolder,'revfig.pdf'));
end

%% Reviewer Figure 1c ISI Histogram, POm
[hist_ctrl,hist_edges,frac_below_1s_ctrl] = collect_isi(db,'Ctrl','spont_test','complex',params);
[hist_cpz,hist_edges,frac_below_1s_cpz] = collect_isi(db,'Cpz','spont_test','complex',params);
figure
hold on
bar(hist_ctrl,'FaceColor',[1 1 1])
bar(hist_cpz,'FaceColor',[0 0 0],'FaceAlpha',0.5)
set(gca,'XTick',(1:3:length(hist_ctrl))+0.5)
set(gca,'Xticklabel',hist_edges(2:3:end))
xlim([1.4 length(hist_ctrl)+1])
xlabel('ISI limit (s)')
ylabel('Fraction of ISIs')
title('POm','FontWeight','normal','Color',params.clr_pom)
ylim([0 0.15])
%legend('Ctrl','Cpz','Location','NorthWest')

%legend boxoff
set(gca,'units','centimeters','position',[3 2 4 3.5]);
if savefigs
     exportgraphics(gcf,fullfile(getdesktopfolder,'revfig.pdf'));
end

%% Reviewer Figure 1d - POm Fraction ISI below 1s Ctrl vs Cpz 
tbl = nj_collect_values( db, conditions,{'frac_below_1s','bursting'},['period_type=spont_test,area=*complex*' add_crit]);
tbl(tbl.frac_below_1s==0,:) = [];
nj_all_data_comparison(tbl,'frac_below_1s','Fraction of ISI below 1s'); %#ok<*UNRCH>
title('POm','FontWeight','normal','Color',params.clr_pom);
if savefigs
     exportgraphics(gcf,fullfile(getdesktopfolder,'revfig.pdf'));
end





%% Helper functions
function [ hist_frac,hist_edges,frac_below_1s] = collect_isi(db,condition,test,area,params)
% Collect ISI from measures

db = db(find_record(db,['condition=' condition]));
hist_edges = params.isi_hist_edges;
mask_below_1s = (hist_edges<=1);
hist_frac = zeros(1,size(hist_edges,2)-1);
frac_below_1s = [];
count = 0;
for i = 1:length(db)
    measures = db(i).measures;
    for j=1:length(measures)
        if isempty(measures(j).isi_hist_frac)
            continue
        end
        if ~strcmp(measures(j).period_type,test)
            continue
        end
        if ~contains(measures(j).area,area)
            continue
        end
        if ~measures(j).good_bc
            continue
        end
        if measures(j).isi_hist_frac(end)>0
            continue % remove silent units
        end

        hist_frac = hist_frac + measures(j).isi_hist_frac;
        count = count+1;

        frac_below_1s(count) = sum(measures(j).isi_hist_frac(mask_below_1s));
    end
end
hist_frac = hist_frac/count;
frac_below_1s_ctrl = frac_below_1s;
end