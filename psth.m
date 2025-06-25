window = [-0.3 0.3]; % look at spike times from 0.3 sec before each event to 1 sec after

myKsDir = uigetdir('\\vs03.herseninstituut.knaw.nl\VS03-AXS-1\NIN212104_Jamann\in_vivo\Neuropixels\SpikeGLX_DATA', 'Select KS directory');
sp = loadKSdir(myKsDir);

% if your events come in different types, like different orientations of a
% visual stimulus, then you can provide those values as "trial groups",
% which will be used to construct a tuning curve. Here we just give a
% vector of all ones. 
eventTimes = pulse_times_all_recordings{2};
trialGroups = ones(size(eventTimes));

psthViewer(sp.st, sp.clu, eventTimes, window, trialGroups);

%%
cluster_num = 492;
spikeTimes = sAP.sCluster(cluster_num).SpikeTimes;
psthViewer(spikeTimes, 10, eventTimes, window, trialGroups);

