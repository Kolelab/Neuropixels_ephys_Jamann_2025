function [start_rec,end_rec,rec_time] = getspontrec(fileList,condition)

%function to find the start and end times of the spontaneous activity
%recordings (between recordings with stimuli or at the start (rec = 0))
% input 1 fileList: list of all files in directory to analyse, has to match
% the vecRecording 
% input 2 condition: Cpz or Ctrl

%preallocate
num_files = numel(fileList);
start_rec = cell(1,num_files);
end_rec = cell(1,num_files);
rec_time = cell(1,num_files);

%% manually make list of all opto recordings
%since its not always directly at the beginning of the recording in my case
% examples of numbers:
% if 0: is at the begginning of recording (untill first stim protocol)
% if 1: during stimulus protocol 1, but stimulus (LED) was switched off
% if 1.5: is between stimulus protocol 1 and 2 

if contains(condition,'Ctrl')
%vecRecording = [0, 7.5, 0, 0, 0, 0, 0, 0, 0, 1.5, 1.5, 0]'; %Ctrl
vecRecording = [0, 7.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]'; %Ctrl, only short recs
else
%vecRecording = [0, 0, 0, 0, 0, 0, 1, 1, 2, 0, -1, 0, 0]'; %Cpz
vecRecording = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2]'; %Cpz, only short recs
end

T_recs = struct2table(fileList);
T_recs.recnumber = vecRecording;


%% Loop for all files to find spont recording activity

for file_idx = 1:num_files
    file_info = fileList(file_idx);
    load(file_info.name)
    
idx_spont_rec = T_recs.recnumber(file_idx);

%find time of spontaneous recording

if idx_spont_rec == -1   %if no stimuli were recorded --> entire file is spontaneous

start_rec{file_idx} = 0; 
end_rec{file_idx} = str2double(sAP.sSources.sMetaNI.fileTimeSecs);   
rec_time{file_idx} = end_rec{file_idx} - start_rec{file_idx};

elseif idx_spont_rec == 0
    
start_rec{file_idx} = 0;     %in case of spontaneous recording before first stimulus
end_rec{file_idx} = sAP.cellBlock{1, 1}.vecStimOnTime(1);  
rec_time{file_idx} = end_rec{file_idx} - start_rec{file_idx};

elseif mod(idx_spont_rec, 1) ~= 0    %if recordings were in between two stim recordings
    
  start_idx = idx_spont_rec-0.5;
  end_idx= start_idx + 1;
  start_rec{file_idx} = sAP.cellBlock{1, start_idx}.vecStimOffTime(end)+26;
  end_rec{file_idx} = sAP.cellBlock{1, end_idx}.vecStimOnTime(1);
  rec_time{file_idx} = end_rec{file_idx} - start_rec{file_idx};

else  %not zero, not integer
      
start_rec{file_idx} = sAP.cellBlock{1, idx_spont_rec}.vecStimOnTime(1) - sAP.cellBlock{1, idx_spont_rec}.dblPrePostWait - sAP.cellBlock{1, idx_spont_rec}.dblPulseWait;
end_rec{file_idx} = sAP.cellBlock{1, idx_spont_rec}.vecStimOffTime(end);
rec_time{file_idx} = end_rec{file_idx} - start_rec{file_idx};

end
end
end

