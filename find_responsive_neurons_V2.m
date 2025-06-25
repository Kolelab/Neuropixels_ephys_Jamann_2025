
%% With this pipeline I want to go through sAP files and look for responsive neurons (whisker and opto) 
% and store this information in sCluster for future analysis

% LOAD FILE and stimulus info

cd(uigetdir);
fileList = dir('*.mat');
num_files = numel(fileList);

file_idx = 13; 
file_info = fileList(file_idx);

full_name = file_info.name;
rec_name = full_name(1:(end-4));


%% load sAP, get stimulus times
[sAP,pulse_times_all_recordings] = getpulsetimes(rec_name);

% Define clusters
% get cluster info
cluster_number = numel([sAP.sCluster.KilosortGood]);  %total clusters in sAP
good_clusters = find(strcmp({sAP.sCluster.bc_unitType}, 'GOOD') == 1);
area_list = cellstr({sAP.sCluster.Area});


%% Find opto responsive

%input opto rec
user_input = inputdlg(strcat('indicate opto test rec for  ',rec_name));
opto_rec = str2num(user_input{1});


for i = 1:cluster_number
zeta(i) = sAP.sCluster(i).ZetaP(opto_rec);
end

resp_clust_opto = find(zeta < 0.05); 
resp_clust_good_opto = intersect(resp_clust_opto,good_clusters);
number_low_zeta_clust = numel(resp_clust_good_opto);

% make raster plot for neuron of choice to visualise response, 
%assign identity 

for Intrespclust_opto = 1:number_low_zeta_clust

idx_cluster = resp_clust_good_opto(Intrespclust_opto);    
vecEventStarts = cell2mat(pulse_times_all_recordings(opto_rec));
vecSpikeTimes = sAP.sCluster(idx_cluster).SpikeTimes;
% f1 = figure;
plotRaster(vecSpikeTimes,vecEventStarts(:,1),0.05,10000)   %Jorrits plotversion
% f1 = (Intrespclust_opto);
% matspikes = getspikes(vecSpikeTimes,vecEventStarts);
% plotSpikeRaster(matspikes, 'PlotType','scatter','XLimForCell',[-0.02 0.1]);
ans_user = inputdlg('unresponsive = 1, responsive = 2, unsure = 3, artefact = 4, art+unresp = 5, art+resp = 6, use-on-artefact = 7, use-off-artefact = 8'); %input responsiveness
dbl_putative_opto = str2double(ans_user{1});
close

sAP.sCluster(idx_cluster).OptoResp = dbl_putative_opto;

end

%% Find whisker responsive

%input whisker rec
user_input = inputdlg(strcat('indicate whisker test rec for  ',rec_name));
whisker_rec = str2num(user_input{1});

for i = 1:cluster_number
whisker_zeta(i) = sAP.sCluster(i).ZetaP(whisker_rec);
end

resp_clust_whisk = find(whisker_zeta < 0.05);  % find putative responsive neurons based on low zeta, for now 0.05 cutoff
resp_clust_good_whisk = intersect(resp_clust_whisk,good_clusters);
num_low_zeta_whisk= numel(resp_clust_good_whisk);
% make raster plot for neuron of choice to visualise response, 
%assign identity 

for Intrespclust_whisk = 1:num_low_zeta_whisk

idx_cluster = resp_clust_good_whisk(Intrespclust_whisk);    
vecEventStarts = cell2mat(pulse_times_all_recordings(whisker_rec));
vecSpikeTimes = sAP.sCluster(idx_cluster).SpikeTimes;
f1 = figure;
plotRaster(vecSpikeTimes,vecEventStarts(:,1),0.05,10000)   %Jorrits plotversion

% f2 = figure(Intrespclust_whisk);
% matspikes = getspikes(vecSpikeTimes,vecEventStarts);
% plotSpikeRaster(matspikes, 'PlotType','scatter','XLimForCell',[-0.02 0.1]);
ans_user = inputdlg('unresponsive = 1, responsive = 2, unsure = 3, artefact = 4, art+unresp = 5, art+resp = 6'); %input responsiveness
dbl_putative_whisker= str2double(ans_user{1});
close

sAP.sCluster(idx_cluster).WhiskerResp = dbl_putative_whisker;

end

%% SAVE

filename = strcat(rec_name , '_resp.mat');
save_loc = uigetdir;
cd(save_loc);
save(filename,'sAP','opto_rec', 'whisker_rec');


