function [SpikeTimes,StimTimes] = find_whiskerstim(files,toggle_whisker_only,toggle_whisker,stimduration)

%FIND whisker stim times from files (either whisker only trials or
%opto/whisker stim combined --> then need to find whisker only trials
%output are two cell arrays with the spike and stimtimes for each files
%(sorted by file)


filenumber = numel(files);
SpikeTimes = cell(1,filenumber);
StimTimes = cell(1,filenumber);

for file = 1:filenumber


%if data from whisker only protocol

if toggle_whisker_only == '1'  
    
%load datatable from excel file
T = readtable(files(file).name);

%get stim times
vec_StimTimes_all = (T.Latency_s_(~isnan(T.Latency_s_)) + (T.Episode_(~isnan(T.Episode_))-1)*20); 

time_diff = diff(vec_StimTimes_all);
idx_one_sec = find(time_diff > 1);

% check if stimuli at 1Hz or different frequencies
%if vec_frequencies(file) == 0
vec_StimTimes = vec_StimTimes_all;   %if only 1Hz stim
%elseif vec_frequencies(file) == 1
% % frequencies: find 1Hz stim only
% time_diff = round(diff(vec_StimTimes_all),0);
% make_string = num2str(time_diff)';
% consecutiveOnesStart = strfind((make_string), '11111');
% consecutiveOnesIndices = [];
% for i = 1:length(consecutiveOnesStart)
%     consecutiveOnesIndices = [consecutiveOnesIndices, consecutiveOnesStart(i):consecutiveOnesStart(i)+4];
% end
% vec_StimTimes = (vec_StimTimes_all(consecutiveOnesIndices)); 
% else
% time_diff = round(diff(vec_StimTimes_all),0);
% make_string = num2str(time_diff)';
% consecutiveOnesStart = strfind((make_string), '1111');
% consecutiveOnesIndices = [];
% for i = 1:length(consecutiveOnesStart)
%     consecutiveOnesIndices = [consecutiveOnesIndices, (consecutiveOnesStart(i)-1):consecutiveOnesStart(i)+3];   
% end
% vec_StimTimes = (vec_StimTimes_all(consecutiveOnesIndices)); 
% end
numStim = numel(vec_StimTimes);

% get spike times
try
vec_SpikeTimes = ((T.Latency_s__1(~isnan(T.Latency_s__1)) + (T.Episode__1(~isnan(T.Episode__1))-1)*20));   %spike times in s
catch
vec_SpikeTimes = [];   % if no spikes were fired
end 
vec_StimEnd = vec_StimTimes + stimduration;


%if  opto + whisker data
else
 
vec_StimTimes_stim1= [];
vec_StimTimes_stim2 = [];
vec_SpikeTimes = [];
    
%load datatable and make lists with stim and spiketimes
T = readtable(files(file).name);     %generate table with data of excel file

vec_StimTimes_stim1 = (T.Latency_s_(~isnan(T.Latency_s_)) + (T.Episode_(~isnan(T.Episode_))-1)*20);  %first stim times (whisker or opto)
vec_StimTimes_stim2 = (T.Latency_s__1(~isnan(T.Latency_s__1)) + (T.Episode__1(~isnan(T.Episode__1))-1)*20);  %second stim times (whisker or opto)
vec_SpikeTimes = (T.Latency_s__2(~isnan(T.Latency_s__2)) + (T.Episode__2(~isnan(T.Episode__2))-1)*20);
if toggle_whisker == '1'
vec_whiskerstim_all = vec_StimTimes_stim1;
else
vec_whiskerstim_all =  vec_StimTimes_stim2;
end


%find the delays between opto and whisker and group them
numStim = numel(vec_StimTimes_stim2);
difference = nan(numStim,1);
for stim2 = 1:numStim
    difference(stim2) = round(vec_StimTimes_stim2(stim2) - vec_StimTimes_stim1(stim2),3);   %delays between stim2 and stim1
end
rounded_difference = round(difference / 0.005) * 0.005;   %round so that it ends in steps of 0.005
[group_idx,Group] = findgroups(rounded_difference);   %split delays in different groups 
numgroups = numel(Group);
whisker_only_trials = find(group_idx ==(find(Group == 1)));
vec_StimTimes = vec_whiskerstim_all(whisker_only_trials);

%vec_SpikeTimes_s = vec_SpikeTimes/1000;
%vec_StimTimes_s = vec_StimTimes/1000;


end

SpikeTimes{file} = vec_SpikeTimes;
StimTimes{file} = vec_StimTimes;

filename_spikes = strcat(files(file).name(1:30), 'spiketimes','.mat');
filename_stims = strcat(files(file).name(1:30), 'stimtimes','.mat');

save(filename_spikes,'vec_SpikeTimes');
save(filename_stims,'vec_StimTimes');

clear T

end











end

