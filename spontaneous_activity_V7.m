%% Script to calculate spontaneous activity from different recorded regions

% This script works in conjuction with the preprocessing pipeline "Acquipix"
% by Jorrit Montijn
% https://github.com/JorritMontijn/Acquipix
% It loads sAP filed (the output of acquipix) and identifies periods of
% spontanoeus activity (user defined). Then it identifies single units
% based on the automated curation tool bombcell 
% (https://github.com/Julie-Fabre/bombcell) and extracts their firing as
% well as bursting behaviour


%%
%first go to directory with AP files (separated Ctrl and Cpz)

cd(uigetdir);
fileList = dir('*.mat');
num_files = numel(fileList);

%% Loop for all files to find start and end of spont recording activity
% edit getspont rec vector according to file list

condition = 'Ctrl';
[start_rec,end_rec,rec_time] = getspontrec(fileList,condition);

%% Find average frequencies of all clusters within spont rec

for file_idx = 1:num_files
    file_info = fileList(file_idx);
    load(file_info.name)
    
try
SampRateImec = sAP.cellBlock{1, 1}.SampRateIM;
catch
SampRateImec = [30000.1167360000];  %for recordings where SampRate missing
end

% get waveform info
numclust = numel([sAP.sCluster.KilosortGood]);  %total clusters in sAP


%% define good clusters, based on bombcell quality metrics 

good_clusters = find(strcmp({sAP.sCluster.bc_unitType}, 'GOOD') == 1); % find good units
area_list = cellstr({sAP.sCluster.Area});

% preallocate
spont_frequency = zeros(1,numel(good_clusters));
vec_Instfreq_spont_bursts = zeros(1,numel(good_clusters));
num_spont_bursts = zeros(1,numel(good_clusters));
vec_burstpermin = zeros(1,numel(good_clusters));

% calculate average spontaneous frequencies  and bursts
for a = 1:numel(good_clusters)
    
    vecspikes = sAP.sCluster(good_clusters(a)).SpikeTimes;  %all spiketimes from one cluster
    vecspikes_spont = vecspikes(find(vecspikes < end_rec{file_idx} & vecspikes > start_rec{file_idx})); %spikes during spontaneous recording
    %average frequency
    num_spikes = numel(vecspikes_spont);   %number of spikes during spont recording
    spont_frequency(a) = num_spikes/rec_time{file_idx};           % frequency (spikes/s)  
    %bursting
    

    if ~isempty(vecspikes_spont)
    spike_threshold = 100; % Minimum spike frequency to be considered a burst (in Hz)
    % Calculate spike intervals
    spike_intervals = diff(vecspikes_spont);
    spike_frequencies = 1 ./ spike_intervals;
    % Find the indices where the spike frequency exceeds the threshold
    burst_indices = find(spike_frequencies > spike_threshold);
    differences = diff(burst_indices);
    burst_remove_indices = (burst_indices(differences == 1))+1; %remove values that belong to a burst with several APs
    burst_indices_removed = burst_indices(~ismember(burst_indices,burst_remove_indices));
    Instfreq_spont_bursts = spike_frequencies(burst_indices_removed); 
    vec_Instfreq_spont_bursts(a) = mean(Instfreq_spont_bursts);
    num_spont_bursts(a) = numel(Instfreq_spont_bursts); 
    vec_burstpermin(a) =  numel(Instfreq_spont_bursts)/(rec_time{file_idx}/60); %bursts per minute, rec time is s

    else
    end

end

spont_frequency_nonzero = spont_frequency(spont_frequency>0);    %find clusters that were actually spiking
idx_nonzero = find(spont_frequency>0);
goodclusters_spiking_idx = good_clusters(idx_nonzero);   %sAP index
number_of_clusters = numel(goodclusters_spiking_idx);

%find bursting cells
cutoff = 10; %define cutoff for number of bursts fired to be counted as bursting unit
vec_Instfreq_spont_bursts_bursting = vec_Instfreq_spont_bursts(num_spont_bursts>cutoff); %if at least X bursts fired
vec_burstpermin_bursting = vec_burstpermin(num_spont_bursts>cutoff);
vec_idx_bursting = find(num_spont_bursts>cutoff);
good_clusters_bursting = good_clusters(vec_idx_bursting);
percentage_bursting = numel(good_clusters_bursting)/numel(good_clusters);

%Make table with mean frequencies per area

areas_animal = area_list(goodclusters_spiking_idx);

% split frequencies of clusters of one animal per area
[AreaIdxList_animal,area_names_animal] = findgroups(area_list(goodclusters_spiking_idx));
[matFrequencies_area_animal] = frequenciesperarea(area_names_animal,AreaIdxList_animal,spont_frequency_nonzero);

%do the same for bursts
%bursting cells
[AreaIdxList_bursting,area_names_bursting] = findgroups(area_list(good_clusters_bursting));
[matbursting_area_animal] = frequenciesperarea(area_names_bursting,AreaIdxList_bursting,vec_Instfreq_spont_bursts_bursting);
[matburstpermin_area_animal] = frequenciesperarea(area_names_bursting,AreaIdxList_bursting,vec_burstpermin_bursting);


% %write excel with frequencies per area per animal
filename = strcat(condition,'_all_clusters_per_animal');  
sheet1 = file_idx;
xlRange1 = 'A1';
xlRange2 = 'B1';
xlRange3 = 'B2';
xlswrite(filename,{file_info.name},sheet1,xlRange1);  
xlswrite(filename,area_names_animal,sheet1,xlRange2);  
xlswrite(filename,matFrequencies_area_animal,sheet1,xlRange3);  

%write excel with bursting per area per animal
filenametwo = strcat(condition,'_bursting_per_animal_min10_100Hz');
sheet1 = file_idx;
xlRange1 = 'A1';
xlRange2 = 'B1';
xlRange3 = 'B2';
xlRange4 = 'A2';
xlswrite(filenametwo,{file_info.name},sheet1,xlRange1);  
xlswrite(filenametwo,cutoff,sheet1,xlRange4); 
xlswrite(filenametwo,area_names_bursting,sheet1,xlRange2);  
xlswrite(filenametwo,matbursting_area_animal,sheet1,xlRange3);  


%write excel with bursting per minute 
filenamethree = strcat(condition,'_bursts_per_minute_minten');
sheet1 = file_idx;
xlRange1 = 'A1';
xlRange2 = 'A2';
xlRange3 = 'A3';
xlRange4 = 'A4';
xlRange5 = 'B1';
xlRange6 = 'B2';

xlswrite(filenamethree,{file_info.name},sheet1,xlRange1);  
xlswrite(filenamethree,cutoff,sheet1,xlRange2); 
xlswrite(filenamethree,{'percentage bursting'},sheet1,xlRange3); 
xlswrite(filenamethree,percentage_bursting,sheet1,xlRange4); 
xlswrite(filenamethree,area_names_bursting,sheet1,xlRange5);  
xlswrite(filenamethree,matburstpermin_area_animal,sheet1,xlRange6);




%concatenate list of all areas/frequencies for all animals
area_array{1,file_idx} = areas_animal;
frequency_array{1,file_idx} = spont_frequency_nonzero;
if file_idx == 1
cat_areas = area_array{1};
cat_frequencies = frequency_array{1};
else
    cat_areas = horzcat(cat_areas,area_array{file_idx});     
    cat_frequencies = horzcat(cat_frequencies,frequency_array{file_idx});
end

end

%% Generate a table with the mean frequencies per area for all animals

[output_path] = uigetdir;
cd (output_path);

[AreaIdxList,area_names] = findgroups(cat_areas);  % find areas in all_areas list
mean_frequency = splitapply(@mean, cat_frequencies, AreaIdxList); % average frequencies per area
SD_frequency = splitapply(@std, cat_frequencies, AreaIdxList); % average SD per area
numel_frequency = splitapply(@numel, cat_frequencies, AreaIdxList); %number of clusters

T_meanfreq_per_area = table(area_names',mean_frequency',SD_frequency',numel_frequency','VariableNames',{'Area','Frequency','SD','n'});
writetable(T_meanfreq_per_area,strcat(condition,'_mean_frequencies_short_recs.xls'),'Range','B2:Z90');

%% Divide cluster frequencies per area over all animals

[matFrequencies_area] = frequenciesperarea(area_names,AreaIdxList,cat_frequencies);

%save excel sheet with frequencies per area
filename = strcat(condition,'_all_clusters_short_recs');  %edit
sheet1 = 1;
sheet2 = 2;
xlRange1 = 'A1';
xlRange2 = 'A2';
xlswrite(filename,area_names,sheet1,xlRange1);  
xlswrite(filename,matFrequencies_area,sheet1,xlRange2);  








