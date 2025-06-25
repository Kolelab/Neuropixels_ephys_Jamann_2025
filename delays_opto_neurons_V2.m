%% open table 
conditions = [{'Ctrl'},{'Cpz'}];

for cnd = 2%:numel(conditions)
condition = conditions{cnd};
    
% condition_ui = inputdlg('Ctrl = 1, Cpz = 2');
% if condition_ui{1} == '1'
% condition = 'Ctrl';
% else
% condition = 'Cpz';
% end

fpath = strcat('\\vs03.herseninstituut.knaw.nl\VS03-AXS-1\NIN212104_Jamann\in_vivo\Neuropixels\Analysis\RecordingProcessor\Final data including waveform\',condition,'\Inspected');
cd(fpath);
fileList = dir('*.mat');
num_files = numel(fileList);

%%

Opto = struct();
Whisker = struct();

for file_idx = 1:num_files


file_info = fileList(file_idx);
load(file_info.name)
name_short = file_info.name(1:(end-9));

[~,pulse_times_all_recordings] = getpulsetimes(condition,name_short);
start_opto_rec = sAP.cellBlock{1, opto_rec}.vecStimOnTime(1) - 0.5;
end_opto_rec = sAP.cellBlock{1, opto_rec}.vecStimOffTime(end) + 0.5; 


%% define responsive neurons
%RN is short for responsive neurons
num_clust_total = numel(sAP.sCluster);

%replace all empty values 
for clust_idx = 1:num_clust_total
    if isempty(sAP.sCluster(clust_idx).OptoResp) 
        sAP.sCluster(clust_idx).OptoResp = 9;
    end
    if isempty(sAP.sCluster(clust_idx).WhiskerResp)
        sAP.sCluster(clust_idx).WhiskerResp = 9;
    end
end

%define neurons
vec_rn_opto = find([sAP.sCluster.OptoResp] == 2); %responsive neurons
vec_rn_unsure_opto = find([sAP.sCluster.OptoResp] == 3); %maybe responsive
vec_rn_mixed_opto = find([sAP.sCluster.OptoResp] == 6); %responsive mixed with artefact
vec_rn_opto_mixed_combined = sort(horzcat(vec_rn_opto, vec_rn_mixed_opto));
vec_rn_opto_and_unsure = sort(horzcat(vec_rn_opto, vec_rn_unsure_opto));
try
vec_rn_artefact = find([sAP.sCluster.OptoResp] == 7); %for alignment
spikes_artefact = sAP.sCluster(vec_rn_artefact(1)).SpikeTimes;
opto_onset = spikes_artefact(find(spikes_artefact < end_opto_rec & spikes_artefact > start_opto_rec));
catch
    vec_rn_artefact = find([sAP.sCluster.OptoResp] == 8);
    spikes_artefact = sAP.sCluster(vec_rn_artefact(1)).SpikeTimes;
opto_onset = (spikes_artefact(find(spikes_artefact < end_opto_rec & spikes_artefact > start_opto_rec)))-0.02;
end
whisker_onset = pulse_times_all_recordings{whisker_rec};
vec_rn_whisker = find([sAP.sCluster.WhiskerResp] == 2);
vec_rn_unsure_whisker = find([sAP.sCluster.WhiskerResp] == 3);

% responsive neurons that are either "Ks good" or have no violations
vec_rn_all_opto = sort(horzcat(vec_rn_opto, vec_rn_unsure_opto, vec_rn_mixed_opto));
vec_rn_all_whisker = sort(horzcat(vec_rn_whisker, vec_rn_unsure_whisker));
vec_good_clusters = find([sAP.sCluster.KilosortGood]==1);
vec_violations = find([sAP.sCluster.Violations1ms]<0.2);   %find out how to do quicker!!
vec_violations2 = find([sAP.sCluster.Violations2ms]<0.2); 
vec_violations_all = unique(sort(horzcat(vec_violations, vec_violations2)));
vec_rn_good_opto = vec_good_clusters(ismember(vec_good_clusters,vec_rn_all_opto));
vec_rn_violations_opto = vec_violations_all(ismember(vec_violations_all,vec_rn_all_opto));
vec_rn_both_opto = unique(sort(horzcat(vec_rn_good_opto, vec_rn_violations_opto)));
vec_rn_good_whisker = vec_good_clusters(ismember(vec_good_clusters,vec_rn_all_whisker));     
vec_rn_violations_whisker = vec_violations_all(ismember(vec_violations_all,vec_rn_all_whisker));   
vec_rn_both_whisker = unique(sort(horzcat(vec_rn_good_whisker, vec_rn_violations_whisker)));

num_rn_opto_good = numel(vec_rn_violations_opto); 
num_rn_whisker_good = numel(vec_rn_violations_whisker);

include_neurons_opto = vec_rn_all_opto; 
include_neurons_whisker = vec_rn_all_whisker; 


%% Neurons in an area
%Area_name = 'Primary somatosensory area';
Area_name = 'Posterior complex';

% find S1 neurons that are opto responsive --> must be L5 CT
for clust = 1:numel([sAP.sCluster])
vec_Area(clust) = contains(sAP.sCluster(clust).Area,Area_name); %edit
end
vec_Area_neurons = find(vec_Area == 1);
vec_Area_opto_neurons = include_neurons_opto(ismember(include_neurons_opto,vec_Area_neurons));
vec_Area_whisker_neurons = include_neurons_whisker(ismember(include_neurons_whisker,vec_Area_neurons));
num_optoneurons_area = numel(vec_Area_opto_neurons);
num_whiskerneurons_area = numel(vec_Area_whisker_neurons);

start = numel(Opto);

%accumulate all relevant info in one big structure for all animals

%opto responsive
if ~isempty(vec_Area_opto_neurons)
    
for Area_loop_idx_o = 1:num_optoneurons_area
    
Opto_idx = vec_Area_opto_neurons(Area_loop_idx_o);
if start == 1
    struct_idx = Area_loop_idx_o;
else
struct_idx = start + Area_loop_idx_o;
end
Opto(struct_idx).Animal = sAP.sCluster(Opto_idx).Subject;
Opto(struct_idx).Cluster = Opto_idx;
Opto(struct_idx).Area = sAP.sCluster(Opto_idx).Area;
Opto(struct_idx).FF = sAP.sCluster(Opto_idx).Spont;
Opto(struct_idx).SpikeTimes = sAP.sCluster(Opto_idx).SpikeTimes;
Opto(struct_idx).Waveform = sAP.sCluster(Opto_idx).Waveform;
Opto(struct_idx).TTPTime = sAP.sCluster(Opto_idx).TTPTime;
Opto(struct_idx).TTPRatio = sAP.sCluster(Opto_idx).TTPRatio;
Opto(struct_idx).Opto = sAP.sCluster(Opto_idx).OptoResp;
Opto(struct_idx).Whisker = sAP.sCluster(Opto_idx).WhiskerResp;
Opto(struct_idx).OptoTimes = opto_onset;

end

else
end

start_w = numel(Whisker);

% whisker responsive
if ~isempty(vec_Area_whisker_neurons)
    
for Area_loop_idx_w = 1:num_whiskerneurons_area
    
Whisker_idx = vec_Area_whisker_neurons(Area_loop_idx_w);
if start_w == 1
    struct_idx = Area_loop_idx_w;
else
struct_idx = start_w + Area_loop_idx_w;
end
Whisker(struct_idx).Animal = sAP.sCluster(Whisker_idx).Subject;
Whisker(struct_idx).Cluster = Whisker_idx;
Whisker(struct_idx).Area = sAP.sCluster(Whisker_idx).Area;
Whisker(struct_idx).FF = sAP.sCluster(Whisker_idx).Spont;
Whisker(struct_idx).SpikeTimes = sAP.sCluster(Whisker_idx).SpikeTimes;
Whisker(struct_idx).Waveform = sAP.sCluster(Whisker_idx).Waveform;
Whisker(struct_idx).TTPTime = sAP.sCluster(Whisker_idx).TTPTime;
Whisker(struct_idx).TTPRatio = sAP.sCluster(Whisker_idx).TTPRatio;
Whisker(struct_idx).Opto = sAP.sCluster(Whisker_idx).OptoResp;
Whisker(struct_idx).Whisker = sAP.sCluster(Whisker_idx).WhiskerResp;
Whisker(struct_idx).WhiskerTimes = whisker_onset;

end

else
end


end

%% GENERATE HEATMAP and averaged response

sOptions.handleFig = -1; %dont plot PEP
vecTraceOrWindow = [-0.02 : 0.005: 0.1];   %bins for heatmap, tracing

%% Opto
num_optoneurons = numel(Opto);

num_bins = numel(vecTraceOrWindow)-1;
heat_matrix_o = nan(num_optoneurons ,num_bins);
%%
for Int_neuron_o = 42%1:num_optoneurons %num_optoneurons 
Plot_neuron = Int_neuron_o;
%spike matrix one neuron
vecSpikes = Opto(Plot_neuron).SpikeTimes;
vecOptoOn = Opto(Plot_neuron).OptoTimes;
matspikes = getspikes(vecSpikes,vecOptoOn);
plotSpikeRaster(matspikes, 'PlotType','scatter','XLimForCell',[-0.02 0.1]);
% ans_user = inputdlg('include = 1, dont include = 0, crop = 2'); %input responsiveness
% dbl_todo = str2double(ans_user{1});
% Opto(Plot_neuron).Use = dbl_todo;
% %clf
% % crop out artefact for those neurons that are contaminated
% if Opto(Plot_neuron).Use == 2
% cutrangeOn = 0.001; %crop around both sides of artefact
% cutrangeOff = 0.0008;
% vecOptoOff = Opto(Plot_neuron).OptoTimes + 0.02;
% for stim = 1:numel(vecOptoOn)
% spikesduringOn = find(vecSpikes > (vecOptoOn(stim)-cutrangeOn) & vecSpikes <  (vecOptoOn(stim)+cutrangeOn));
% spikesduringOff = find(vecSpikes > (vecOptoOff(stim)-cutrangeOff) & vecSpikes < (vecOptoOff(stim)+cutrangeOff));
% vecSpikes(spikesduringOn) = [];
% vecSpikes(spikesduringOff) = [];
% end
% matspikes = getspikes(vecSpikes,vecOptoOn);
% plotSpikeRaster(matspikes, 'PlotType','scatter','XLimForCell',[-0.02 0.1]);
% ans_user_2 = inputdlg('now include? yes == 1, no == 9'); %input if properly cropped
% dbl_todo_2 = str2double(ans_user_2{1});
% Opto(Plot_neuron).Use = dbl_todo_2;
% Opto(Plot_neuron).SpikeTimes = vecSpikes;
% clf
% else
% end
end

%%
Opto_use = find([Opto.Use]==1);
num_opto_use = numel(Opto_use);
heat_matrix_o = nan(num_opto_use,num_bins);

for Int_Opto_use = 1:num_opto_use
    Idx_un = Opto_use(Int_Opto_use);
    vecSpikes_un = Opto(Idx_un).SpikeTimes;
    vecOptoOn_un = Opto(Idx_un).OptoTimes;
    num_trials = numel(vecOptoOn_un);
[vecMean,vecSEM,vecWindowBinCenters,matPET] = doPEP(vecSpikes_un,vecTraceOrWindow,vecOptoOn_un,sOptions);

vecz_o = zscore(matPET(:));   %z-score the spike matrix
zmatrix_o = reshape(vecz_o,[num_trials,num_bins]);
vecMean_z = nanmean(zmatrix_o,1);
sd_z = nanstd(zmatrix_o,1);

spike_mean = nanmean(matPET,1); %first mean, then z-score 
mean_zscored = zscore(spike_mean);
heat_matrix_o(Int_Opto_use,:) = mean_zscored; 
mean_vec = nanmean(heat_matrix_o,1);
sd_vec = nanstd(heat_matrix_o,1);
end

[~,Maxidx] = max(heat_matrix_o, [], 2);
[~, sortIdx_o] = sort(Maxidx,'ascend');
sorted_matrix_o  = heat_matrix_o(sortIdx_o, :);


% f = figure(cnd);
subplot(2,2,1);
x = vecTraceOrWindow(1:num_bins);
y = [1:1:(numel(num_optoneurons))];
imagesc(x,y,sorted_matrix_o)
colorbar

subplot(2,2,2);
errorfill(x,mean_vec,sd_vec);

for i = 1:numel(Maxidx)
peaktime(i) = vecTraceOrWindow(Maxidx(i));
sd(i)  = std(sorted_matrix_o(i,:));
end

%% Whisker

num_wn = numel(Whisker);
vecTraceOrWindow_w = [-0.05 : 0.001: 0.1];   %5ms bins
num_bins = numel(vecTraceOrWindow_w)-1;
%heat_matrix_w = nan(num_wn ,num_bins);

for Int_neuron_w = 1:num_wn 
vecSpikes = Whisker(Int_neuron_w).SpikeTimes;
vecWhiskerOn = Whisker(Int_neuron_w).WhiskerTimes;
matspikes = getspikes(vecSpikes,vecWhiskerOn);
plotSpikeRaster(matspikes, 'PlotType','scatter','XLimForCell',[-0.02 0.1]);
ans_user = inputdlg('include = 1, dont include = 0'); %input responsiveness
dbl_todo = str2double(ans_user{1});
Whisker(Int_neuron_w).Use = dbl_todo;
clf
end


%%


Area_name = 'Primary somatosensory area';
%Area_name = 'Posterior complex';
num_wn = numel(Whisker);
vecTraceOrWindow_w = [-0.05 : 0.001: 0.1]; 
num_bins = numel(vecTraceOrWindow_w)-1;

num_clust = numel([Whisker.Opto]);
vec_Area = nan(num_clust,1);
% find S1 neurons that are opto responsive --> must be L5 CT
for clust = 1:num_clust
vec_Area(clust) = contains(Whisker(clust).Area,Area_name); %edit
end

vec_Area_neurons = find(vec_Area == 1);

Whisker_use = find([Whisker.Use]==1);
Whisker_use_L5 = Whisker_use(ismember(Whisker_use,vec_Area_neurons)); 
num_whisker_use = numel(Whisker_use_L5);
heat_matrix_w = nan(num_whisker_use,num_bins);


for Int_whisker_use = 1:num_whisker_use
    Idx_un = Whisker_use_L5(Int_whisker_use);
    vecSpikes_un_w = Whisker(Idx_un).SpikeTimes;
    vecWhiskerOn_un_w = Whisker(Idx_un).WhiskerTimes;
    
    num_trials_w = numel(vecWhiskerOn_un_w);

[vecMean,vecSEM,vecWindowBinCenters,matPET] = doPEP(vecSpikes_un_w,vecTraceOrWindow_w,vecWhiskerOn_un_w,sOptions);

veczW = zscore(matPET(:));   %z-score the spike matrix
zmatrixW = reshape(veczW,[num_trials_w,num_bins]);
vecMean_z = nanmean(zmatrixW ,1);
sd_vec_z = nanstd(zmatrixW,1);

spike_mean = nanmean(matPET,1); %first mean, then z-score
mean_zscored = zscore(spike_mean);
heat_matrix_w(Int_whisker_use,:) = mean_zscored; 
mean_vec = nanmean(heat_matrix_w,1);
sd_vec = nanstd(heat_matrix_w,1);
end

animals = [Whisker(:).Animal];
[M,Maxidx_w] = max(heat_matrix_w, [], 2);
[~, sortIdx_w] = sort(Maxidx_w,'ascend');
sorted_matrix_w  = heat_matrix_w(sortIdx_w, :);

subplot(2,2,3);
x = vecTraceOrWindow_w(1:num_bins);
y = [1:1:num_whisker_use];
imagesc(x,y,sorted_matrix_w)
colorbar

subplot(2,2,4);
errorfill(x,mean_vec,sd_vec);

for i = 1:numel(Maxidx_w)
peaktime(i) = vecTraceOrWindow_w(Maxidx_w(i));
end

% vars = {'Opto','Whisker'};
% clear(vars{:})

end

