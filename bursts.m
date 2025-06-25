%% detect bursts

cd(uigetdir);
fileList = dir('*.mat');
num_files = numel(fileList);

condition = 'Ctrl';
[start_rec,end_rec,rec_time] = getspontrec(fileList,condition);

for file_idx = 1:num_files

file_info = fileList(file_idx);
load(file_info.name)

%% find PoM neurons

index_PoM = find(strcmp({sAP.sCluster.Area}, 'Posterior complex of the thalamus')==1);
index_good = find([sAP.sCluster.Violations1ms] <0.2 | [sAP.sCluster.Violations2ms] <0.2);

index_good_Pom = index_good(ismember(index_good, index_PoM));

%% plot all their activity binned

vec_depth = [sAP.sCluster.Depth];
depths = vec_depth(index_good_Pom);
%select time-window you want to plot (in s)

time_start = start_rec{file_idx};
time_end = end_rec{file_idx};
ISI_cell = {};
spike_mat = nan(1,1);

for PoMidx = 1:numel(index_good_Pom)
    
    vecSpikes = sAP.sCluster(index_good_Pom(PoMidx)).SpikeTimes;
    spikes_in_window = vecSpikes(find(vecSpikes > time_start & vecSpikes < time_end));
    number_spikes = numel(spikes_in_window);
    spike_mat(PoMidx,1:number_spikes) = spikes_in_window;
    if number_spikes > 2
    for spike = 2:number_spikes
    ISI(spike-1) = spikes_in_window(spike) - spikes_in_window(spike-1);
    ISI_cell{PoMidx} = ISI;
    end
    else
    ISI_cell{PoMidx} = nan;
    end
  
end


%% find bursts
vec_idx_bursts = nan;
num_burst = nan;
for i = 1:numel(index_good_Pom)
vec_idx_bursts = find(ISI_cell{i} < 0.01);
bursts{i} = vec_idx_bursts; 
num_burst(i) = numel(vec_idx_bursts)/rec_time{file_idx}; %number of bursts/s
end
%%

%  y2 = depths;
%  x2 = spike_mat;  
%  plot(x2,y2,'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor',[0 0 0],'MarkerSize',2)
%  xlabel('time (s)')
%  xlim([time_start time_end])
% %  ylim([1900 1905])
if ~isnan(num_burst)==1
average_bursts{file_idx} = nanmean(num_burst);
number_bursts{file_idx} = num_burst(:);
else
end

end
 %%
 
save('bursts_Cpz.mat','average_bursts');
