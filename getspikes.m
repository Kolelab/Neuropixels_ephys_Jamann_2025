function [spikes] = getspikes(vecSpikeTimes,vecEventStarts)

%SPIKES  find all spikes of one cluster of predefined opto recording,
%output relative delays to onset of stimulus in cell array called spikes
%input1: idx_cluster = index of cluster 
%input2 sAP (Neuropixel data synthesis)
%input 3: opto_rec = index of recording
% input 4: vecSpikes_opto = spikes of artefact cluster for alignmnet of
% onset

%vecSpikeTimes = sAP.sCluster(idx_cluster).SpikeTimes;  
% start_rec = sAP.cellBlock{1, opto_rec}.vecStimOnTime(1);
% end_rec = sAP.cellBlock{1, opto_rec}.vecStimOnTime(end);
% vecEventStarts = vecSpikes_opto(find(vecSpikes_opto < end_rec & vecSpikes_opto > start_rec));
num_events = numel(vecEventStarts);
spikes = cell(num_events,1);  %all the spikes in the chosen window separated by trial 
plot_window = 0.1;  %window for raster plot (before and after onset)

for idxtrial = 1:num_events
    startTime = vecEventStarts(idxtrial)-plot_window;
    endTime = startTime + (plot_window*2);
    idx_spikes_in_window = find(vecSpikeTimes < endTime & vecSpikeTimes > startTime); %find spikes in the indicated time window around the pulse onset
    aligned_times = vecSpikeTimes(idx_spikes_in_window)-startTime-plot_window; %align to stim onset (t=0)
    if isempty(aligned_times) == 1
         spikes{idxtrial} = nan;   %will be outside plotting window, but maybe find a way to leave empty??
    else
    spikes{idxtrial,1} = aligned_times';
    end

end

