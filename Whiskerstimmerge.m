%This script analyses exported data from axograph saved as .xls files.
%It finds spike and whisker stim times (either from protocols with whisker
%stimulation only or opto+whisker stimulation
%It analyses the delay, variance and probability of spiking in response to
%a series of stimuli, and plots the PSTH so that the user
%can vusally inspect the cells data
%THen it calculates zeta to determine if the cell was significantly
%modulated by the stimulus
%The results are saved in an excel file
%% Find folders

stimduration = 0.020;
addition = 0.050; %indicate the whole window to look for spikes (on top of the stim length) (in s)
spont_period = 0.002; % if earlier than this period after onset then spontaneous

selectedFolder = uigetdir('Select a folder');

% List all subfolders
subfolders = dir(fullfile(selectedFolder, '**', ''));
subfolders = subfolders([subfolders.isdir] & ~ismember({subfolders.name}, {'.', '..'}));

%% Load spike and stimtimes from the different folders and concatenate them (with a delay)

% Loop through subfolders
for cellname = 1:length(subfolders)
    currentSubfolder = fullfile(subfolders(cellname).folder, subfolders(cellname).name);
    addpath(currentSubfolder);

    % List all files in the current subfolder
    filesInSubfolder = dir(fullfile(currentSubfolder, '*.mat'));
    numfiles = numel(filesInSubfolder); 

    spiketimes_first = load(filesInSubfolder(1).name);
    stimtimes_first = load(filesInSubfolder(2).name);
    spiketimes_second = load(filesInSubfolder(3).name);
    stimtimes_second = load(filesInSubfolder(4).name);
    last_time_first_protocol = max(spiketimes_first.vec_SpikeTimes(end),stimtimes_first.vec_StimTimes(end));
    spiketimes_second_new = spiketimes_second.vec_SpikeTimes + last_time_first_protocol + 2;
    stimtimes_second_new = stimtimes_second.vec_StimTimes + + last_time_first_protocol + 2;
    last_time_second_protocol_new = max(spiketimes_second_new(end),stimtimes_second_new(end));
    if numfiles == 6
    spiketimes_third = load(filesInSubfolder(5).name);
    stimtimes_third = load(filesInSubfolder(6).name);  
    spiketimes_third_new = spiketimes_third.vec_SpikeTimes + last_time_second_protocol_new + 2;
    stimtimes_third_new = stimtimes_third.vec_StimTimes + last_time_second_protocol_new + 2;
    else
    spiketimes_third_new = [];
    stimtimes_third_new = [];
    end
    
    SpikeTimes{cellname} = cat(1,spiketimes_first.vec_SpikeTimes,spiketimes_second_new,spiketimes_third_new);
    StimTimes{cellname}  = cat(1,stimtimes_first.vec_StimTimes,stimtimes_second_new,stimtimes_third_new);

end


%% Extract delay, variance, probability, run psth and zeta

excelsheet_filename = '240123_merged_delays';
filenumber = numel(subfolders);

ALL_first_Delays_cat = nan(1000,filenumber);
ALL_Delays_cat = nan(1000,filenumber);

for file = 1:filenumber

vec_SpikeTimes = cell2mat(SpikeTimes(1,file));  
vec_StimTimes = cell2mat(StimTimes(1,file));  
    
stimduration = 0.020; %indicate length of opto stim (in s)
addition = 0.050; %indicate the whole window to look for spikes (on top of the stim length) (in s)
spont_period = 0.002; % if earlier than this period after onset then spontaneous

[probability,first_delay,first_variance,mean_number_of_spikes,num_trials,mean_inst_freq,number_of_bursts,ALL_Delays,ALL_first_Delays] = extract_delay_variance_probability(vec_SpikeTimes,vec_StimTimes,stimduration,addition,spont_period);

%Cumulate average values for each cell(file) in new vector
vec_probability(file) = probability;
vec_mean_delay(file) = first_delay;
vec_variance_delay(file)= first_variance;
vec_mean_number_of_spikes(file) = mean_number_of_spikes;
vec_num_trials(file) = num_trials;
headers(file) = cellstr(subfolders(file).name);
vec_mean_inst_freq(file) = mean_inst_freq;  
vec_number_of_bursts(file) = number_of_bursts;
ALL_Delays_cat(1:numel(ALL_Delays),file) = ALL_Delays';
ALL_first_Delays_cat(1:numel(ALL_first_Delays),file) = ALL_first_Delays';




%% PLot PSTH and decide if responsive

   
% psth = mpsth(vec_SpikeTimes,vec_StimTimes,'fr',0,'binsz',5,'pre',200,'post',500,'chart',1);
% 
% filename_psth= strcat((subfolders(file).name),'_psth');
% savefig(filename_psth)
% answer = inputdlg('none = 0, early = 1, late = 2 depressed = 3, unsure = 4');
% vec_response_inspected(file) = answer; 
% close


%% zetatest
% now with preceding baseline

% dblBaselineDuration = 0.1;
% dblUseMaxDur = 0.5;
% 
% 
% [Zeta,ZetaP,Zetapeak] = zeta_baseline(dblBaselineDuration,dblUseMaxDur,vec_StimTimes,vec_SpikeTimes);
% 
% vec_Zeta(1,file) = Zeta;
% vec_ZetaP(1,file) = ZetaP;
% vec_Zetapeak(1,file) = Zetapeak;
% 
% filename_zeta = strcat((subfolders(file).name),'_zeta');
% savefig(filename_zeta)
% 
% close

end

%% SAVE EXCEL SHEET
%Save output to excel sheet (created in the beginning)

tdir=uigetdir;  %saving location
cd(tdir)

text1 = cellstr('mean delay');
text2 = cellstr('variance');
text3 = cellstr('mean number of spikes');
text4 = cellstr('probability');
text5 = cellstr('number of trials');
text6 = cellstr('burst frequency');
text7 = cellstr('burst number');
text8 = cellstr('responsive');
text9 = cellstr('Zeta');
text10 = cellstr('Zeta p value');
text11 = cellstr('Zeta peak delay');
text12 = cellstr('ALL delays');
text13 = cellstr('ALL first delays');


xlswrite(excelsheet_filename,text1,1,'A2');
xlswrite(excelsheet_filename,text2,1,'A3');
xlswrite(excelsheet_filename,text3,1,'A4');
xlswrite(excelsheet_filename,text4,1,'A5');
xlswrite(excelsheet_filename,text5,1,'A6');
xlswrite(excelsheet_filename,text6,1,'A7');
xlswrite(excelsheet_filename,text7,1,'A8');
xlswrite(excelsheet_filename,text8,1,'A9');
xlswrite(excelsheet_filename,text9,1,'A10');
xlswrite(excelsheet_filename,text10,1,'A11');
xlswrite(excelsheet_filename,text11,1,'A12');
xlswrite(excelsheet_filename,headers,1,'B1');
xlswrite(excelsheet_filename,vec_mean_delay,1,'B2');
xlswrite(excelsheet_filename,vec_variance_delay,1,'B3');
xlswrite(excelsheet_filename,vec_mean_number_of_spikes,1,'B4');
xlswrite(excelsheet_filename,vec_probability,1,'B5');
xlswrite(excelsheet_filename,vec_num_trials,1,'B6');
xlswrite(excelsheet_filename,vec_mean_inst_freq,1,'B7');
xlswrite(excelsheet_filename,vec_number_of_bursts,1,'B8');
%xlswrite(excelsheet_filename,vec_response_inspected,1,'B9');
% xlswrite(excelsheet_filename,vec_Zeta,1,'B10');
% xlswrite(excelsheet_filename,vec_ZetaP,1,'B11');
% xlswrite(excelsheet_filename,vec_Zetapeak,1,'B12');
xlswrite(excelsheet_filename,text12,2,'A2');
xlswrite(excelsheet_filename,text13,3,'A2');
xlswrite(excelsheet_filename,headers,2,'B1');
xlswrite(excelsheet_filename,headers,3,'B1');
xlswrite(excelsheet_filename,ALL_Delays_cat,2,'B2');
xlswrite(excelsheet_filename,ALL_first_Delays_cat,3,'B2');

