%% Script to plot spontaneous acitivy from different recorded regions
%first go to directory with AP files (separated Ctrl and Cpz)

cd(uigetdir);
fileList = dir('*.mat');
num_files = numel(fileList);

condition = 'Cpz';
[start_rec,end_rec,rec_time] = getspontrec(fileList,condition);

%% select recording you want to plot, load sAP

for file_idx = 1:num_files
file_info = fileList(file_idx);
load(file_info.name)

%% Find areas of this recording 

good_clusters = find([sAP.sCluster.KilosortGood]==1);
area_list = cellstr({sAP.sCluster.Area});  %cell array of all areas in sAP
depth_list = cell2mat({sAP.sCluster.Depth}); % cell array of all depths in sAP
[A,areas] = findgroups(area_list); %which areas are recorded

%% Find all good neurons in certain areas of interest (later can be combined in loop)

%PoM
PoMstr = 'Posterior complex';
PoM_idx = find(contains(areas,PoMstr)==1);
if ~isempty(PoM_idx)
PoM_all_cells = find(A == PoM_idx);
PoM_good_cells = PoM_all_cells(find(ismember(PoM_all_cells,good_clusters)==1));
else
end

%S1
S1str = 'Primary somatosensory';
S1_idx = find(contains(areas,S1str)==1);
S1_all_cells = cell(1,numel(S1_idx));
S1_all_cat = [];
if ~isempty(S1_idx)
for layer_idx = 1:numel(S1_idx)
S1_all_cells{layer_idx} = find(A == S1_idx(layer_idx));
end
S1_cat = horzcat(cell2mat(S1_all_cells));
S1_all_cells = S1_cat(find(ismember(S1_cat,good_clusters)==1));
else
end

%CA1
CA1str = 'CA1';
CA1_idx = find(contains(areas,CA1str)==1);
if ~isempty(CA1_idx)
CA1_all_cells = find(A == CA1_idx);
CA1_good_cells = CA1_all_cells(find(ismember(CA1_all_cells,good_clusters)==1));
else disp 'no CA1 cells'
end

%CA2
CA2str = 'CA2';
CA2_idx = find(contains(areas,CA2str)==1);
if ~isempty(CA2_idx)
CA2_all_cells = find(A == CA2_idx);
CA2_good_cells = CA2_all_cells(find(ismember(CA2_all_cells,good_clusters)==1));
else disp 'no CA2 cells'
end

%CA3
CA3str = 'CA3';
CA3_idx = find(contains(areas,CA3str)==1);
if ~isempty(CA3_idx)
CA3_all_cells = find(A == CA3_idx);
CA3_good_cells = CA3_all_cells(find(ismember(CA3_all_cells,good_clusters)==1));
else disp 'no CA3 cells'
end

%LP
LPstr = 'Lateral posterior';
LP_idx = find(contains(areas,LPstr)==1);
if ~isempty(LP_idx)
LP_all_cells = find(A == LP_idx);
LP_good_cells = LP_all_cells(find(ismember(LP_all_cells,good_clusters)==1));
else disp 'no LP cells'
end

%% Plot spikes per area of interest during spontaneous activity rec

start_rec_offset = cell2mat(start_rec(file_idx)) + 0; % indicate start of plot window
end_window =  cell2mat(end_rec(file_idx));  %if you want to plot the entire recording
%end_window = start_rec_offset + 150;   % if only want to plot a part

if exist('PoM_good_cells') == 1
spikeTimes_spont_PoM = cell(1,numel(PoM_good_cells));
for PoM_cell_idx = 1:numel(PoM_good_cells)
    
    spikeTimes = sAP.sCluster(PoM_good_cells(PoM_cell_idx)).SpikeTimes;
    spikeTimes_spont_PoM{PoM_cell_idx} = spikeTimes(find(spikeTimes < end_window & spikeTimes > start_rec_offset));
end
else
end

if exist('S1_all_cells') == 1
spikeTimes_spont_S1 = cell(1,numel(S1_all_cells));
for S1_cell_idx = 1:numel(S1_all_cells)
    spikeTimes = sAP.sCluster(S1_all_cells(S1_cell_idx)).SpikeTimes;
    spikeTimes_spont_S1{S1_cell_idx} = spikeTimes(find(spikeTimes < end_window & spikeTimes > start_rec_offset));
end
else
end

if exist('CA3_all_cells') == 1
spikeTimes_spont_CA3 = cell(1,numel(CA3_all_cells));
for CA3_cell_idx = 1:numel(CA3_all_cells)

    spikeTimes = sAP.sCluster(CA3_all_cells(CA3_cell_idx)).SpikeTimes;
    spikeTimes_spont_CA3{CA3_cell_idx} = spikeTimes(find(spikeTimes < end_window & spikeTimes > start_rec_offset));
end
else 
end

if exist('LP_all_cells') == 1
spikeTimes_spont_LP = cell(1,numel(LP_all_cells));
for LP_cell_idx = 1:numel(LP_all_cells)

    spikeTimes = sAP.sCluster(LP_all_cells(LP_cell_idx)).SpikeTimes;
    spikeTimes_spont_LP{LP_cell_idx} = spikeTimes(find(spikeTimes < end_window & spikeTimes > start_rec_offset));
end
else
end

y1 = (depth_list(PoM_good_cells))*-1;
y2 = (depth_list(S1_all_cells))*-1;
y3 = (depth_list(CA3_all_cells))*-1;
y4 = (depth_list(LP_all_cells))*-1;

for plot_idx = 1:numel(spikeTimes_spont_PoM)
    if ~isempty (spikeTimes_spont_PoM{plot_idx})
plot(spikeTimes_spont_PoM{plot_idx},y1(plot_idx),'r.')
    else
    end
hold on
end
hold on

for plot_idx = 1:numel(spikeTimes_spont_S1)
    if ~isempty (spikeTimes_spont_S1{plot_idx})
plot(spikeTimes_spont_S1{plot_idx},y2(plot_idx),'b.')
hold on
    else disp 'x'
    end
end

hold on

for plot_idx = 1:numel(spikeTimes_spont_CA3)
    if ~isempty (spikeTimes_spont_CA3{plot_idx})
plot(spikeTimes_spont_CA3{plot_idx},y3(plot_idx),'g.')
hold on
    else disp 'x'
    end
end

hold on

for plot_idx = 1:numel(spikeTimes_spont_LP)
    if ~isempty (spikeTimes_spont_LP{plot_idx})
plot(spikeTimes_spont_LP{plot_idx},y4(plot_idx),'b.')
hold on
    else
    end
end

clearvars 'PoM_good_cells' 'S1_all_cells'
%% Calculate binned spikes per area (PoM)

%PoM

if exist('spikeTimes_spont_PoM') == 1
 bin_width = 0.1;
 number_bin = round(rec_time{file_idx}/bin_width);
 num_cells = numel(spikeTimes_spont_PoM);
 num_spikes = zeros(num_cells,number_bin);

 for bin_idx = 1:number_bin
     
     start_bin = start_rec{file_idx} + bin_width*(bin_idx-1);
     end_bin = start_rec{file_idx} + bin_width*bin_idx;
     for cell_idx = 1:numel(spikeTimes_spont_PoM)
         spikes_cell = cell2mat(spikeTimes_spont_PoM(cell_idx));
     num_spikes(cell_idx,bin_idx) = numel(find(spikes_cell < end_bin & spikes_cell > start_bin));
     end   
 end

num_all_spikes{file_idx} = sum(num_spikes,1);
num_cells_all{file_idx} = num_cells;

else
    num_all_spikes{file_idx} = [];

% plot([1:1:numel(num_all_spikes)],num_all_spikes)
% xlim([0 1500]);
end
clearvars 'spikeTimes_spont_PoM'
end

%% FFT

for fft_idx = 1:num_files
    if ~isempty(num_all_spikes{fft_idx})
Fs = 1/bin_width;            % Sampling frequency                    
T = 1/Fs;                    % Sampling period    
L = rec_time{fft_idx};
%L = all_rec_times{fft_idx};       % Length of signal

Y = fft(num_all_spikes{fft_idx});
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
P1_array{fft_idx} = P1;
f{fft_idx} = Fs*(0:(L/2))/L; 
    else f{fft_idx} = [];
        P1_array{fft_idx} = [];
    end 
end

%%
for plot_idx = 1:num_files
    
plot(f{plot_idx},P1_array{plot_idx}) 
xlim([0 5])
ylim([0 20])
xlabel('f (Hz)')
ylabel('|P1(f)|')
hold on
end

%% averaged over recordings

x_vec = [0:0.01:5];
avg_matrix = zeros(num_files,numel(x_vec));

for avg_idx = 1:num_files
    
    for plot_idx2 = 1:(numel(x_vec)-1)
       top = x_vec(plot_idx2+1);  
       bottom = x_vec(plot_idx2);
       f_idx = find(f{avg_idx} < top & f{avg_idx} >= bottom);
       vec_values = P1_array{avg_idx};
       if ~isempty(f_idx)
       avg_matrix(avg_idx,plot_idx2) = vec_values(f_idx(1));
       else
       avg_matrix(avg_idx,plot_idx2) = 0;
       end
    end
end

for plot_idx3 = 1:numel(x_vec)
    
    average(plot_idx3) = mean(nonzeros(avg_matrix(:,plot_idx3)));
    SD(plot_idx3) = std(nonzeros(avg_matrix(:,plot_idx3)));
    n(plot_idx3) = std(nonzeros(avg_matrix(:,plot_idx3)));
    
end


%%
plot(x_vec,average)

[output_path] = uigetdir;
cd (output_path);

T_average = table(x_vec',average',SD',n');
writetable(T_average,'Frequencies_Cpz_PoM.xls','Range','B2')
