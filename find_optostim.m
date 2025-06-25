function [SpikeTimes,StimTimes,totaltime] = find_optostim(files,stimduration)

%Find the opto stim times and spiketimes from excel files 
%output are two cell arrays with the spike and stimtimes for each files
%(sorted by file)


filenumber = numel(files);
SpikeTimes = cell(1,filenumber);
StimTimes = cell(1,filenumber);

for file = 1:filenumber

%fprintf(num2str(file));

%load datatable from excel file
T = readtable(files(file).name);

%get stim times
vec_StimTimes = (T.Latency_s_(~isnan(T.Latency_s_)) + (T.Episode_(~isnan(T.Episode_))-1)*20); 
numStim = numel(vec_StimTimes);

% get spike times
try
vec_SpikeTimes = ((T.Latency_s__1(~isnan(T.Latency_s__1)) + (T.Episode__1(~isnan(T.Episode__1))-1)*20));   %spike times in s
catch
vec_SpikeTimes = [];   % if no spikes were fired
end 
vec_StimEnd = vec_StimTimes + stimduration;
SpikeTimes{file} = vec_SpikeTimes;
StimTimes{file} = vec_StimTimes;


numepisodes = max(T.Episode_); 
totaltime = (numepisodes*20)/ 60;  %in min

end

end












