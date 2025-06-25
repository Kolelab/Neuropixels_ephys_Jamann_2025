function [start_rec,end_rec,rec_time] = getspontrec_singlefile(idx_spont_rec,sAP)

%function to find the start and end times of the sponetaneous activity
%recordings
% input 1 = userinput which spont rec 
%% 
%find time of spontaneous recording

if idx_spont_rec == 0
    
start_rec = 0;     %in case of spontaneous recording before first stimulus
try
end_rec = sAP.cellBlock{1, 1}.vecStimOnTime(1); 
catch
end_rec = 3600;
end
rec_time = end_rec - start_rec;

elseif isinteger(idx_spont_rec) == 0   %if recordings were in between two stim recordings
    
  start_idx = idx_spont_rec-0.5;
  end_idx= start_idx + 1;
  start_rec = sAP.cellBlock{1, start_idx}.vecStimOffTime(end)+26;
  end_rec = sAP.cellBlock{1, end_idx}.vecStimOnTime(1);
  rec_time = end_rec - start_rec;

else
      
start_rec = sAP.cellBlock{1, idx_spont_rec}.vecStimOnTime(1) - sAP.cellBlock{1, idx_spont_rec}.dblPrePostWait - sAP.cellBlock{1, idx_spont_rec}.dblPulseWait;
end_rec = sAP.cellBlock{1, idx_spont_rec}.vecStimOffTime(end);
rec_time = end_rec - start_rec;

end


