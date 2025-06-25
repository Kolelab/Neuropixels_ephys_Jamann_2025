%% This script anlayses the probability during paired opto and whisker stimulation

%% FIND FILES
%find the files with values of : whisker stimulation times, opto
%stimulation times and spiketimes, which were extracted with Axograph from
%the raw data

%User locates target directory with exel files
tdir=uigetdir;
cd(tdir)

files = dir('*.xlsx');

filenumber = numel(files);
stimduration = 0.020; %indicate length of opto stim (s)
addition = 0.15; %indicate the whole window (s) to look for spikes 
addition_all = 0.15; %look for spikes during 100ms window --> both stimuli plus spont spiking!

prompt = 'Whisker first? Yes = 1, No = 0';
toggle_whisker = cell2mat(inputdlg(prompt));
excelsheet_filename = '240516_opto_1st_150mswindow.xlsx';

%% LOAD DATA AND DETECT STIMULUS TIMES

cellname_excel = 'all_cells_num_spikes_opto_first_150ms_240516';

matrix_probabilities_opto = nan(12,filenumber);
matrix_probabilities_whisker = nan(12,filenumber);
matrix_probabilities_all = nan(12,filenumber);
matrix_spikenumber_opto = nan(12,filenumber);
matrix_spikenumber_whisker = nan(12,filenumber);
matrix_spikenumber_all = nan(12,filenumber);
matrix_groups = nan(12,filenumber);

for file = 1:filenumber   %loop for all files (neurons)

vec_StimTimes_stim1= [];
vec_StimTimes_stim2 = [];
vec_SpikeTimes = [];
    
%load datatable and make lists with stim and spiketimes
T = readtable(files(file).name);     %generate table with data of excel file
vec_StimTimes_stim1 = (T.Latency_s_(~isnan(T.Latency_s_)) + (T.Episode_(~isnan(T.Episode_))-1)*20);  %first stim times (whisker or opto)
vec_StimTimes_stim2 = (T.Latency_s__1(~isnan(T.Latency_s__1)) + (T.Episode__1(~isnan(T.Episode__1))-1)*20);  %second stim times (whisker or opto)
vec_SpikeTimes = (T.Latency_s__2(~isnan(T.Latency_s__2)) + (T.Episode__2(~isnan(T.Episode__2))-1)*20);          % spike times

%for all opto and whisker stims, as well as both
% check if spiking (probability) and if yes, 
% how many spikes within a window (number of spikes)

%preallocate
numStim = numel(vec_StimTimes_stim2); %number of stimuli
Spiking_opto = nan(1,numStim);
Num_spikes_opto = nan(1,numStim);
Spiking_whisker = nan(1,numStim);
Num_spikes_whisker = nan(1,numStim);
Spiking_all = nan(1,numStim);
Num_spikes_all = nan(1,numStim);

%find spikes during start of the stimuli (opto/whisker/both="all")

for stim = 1:numStim
    if toggle_whisker == '1'
    vec_SpikesThisTrial_opto = vec_SpikeTimes(vec_SpikeTimes > vec_StimTimes_stim2(stim) & vec_SpikeTimes < (vec_StimTimes_stim2(stim) + addition)); %whisker first
    vec_SpikesThisTrial_whisker = vec_SpikeTimes(vec_SpikeTimes > vec_StimTimes_stim1(stim) & vec_SpikeTimes < (vec_StimTimes_stim1(stim) + addition));
    vec_SpikesThisTrial_all = vec_SpikeTimes(vec_SpikeTimes > vec_StimTimes_stim1(stim) & vec_SpikeTimes < (vec_StimTimes_stim1(stim) + addition_all)); % look during entire window (both wisker and opto)
    else
    vec_SpikesThisTrial_opto = vec_SpikeTimes(vec_SpikeTimes > vec_StimTimes_stim1(stim) & vec_SpikeTimes < (vec_StimTimes_stim1(stim) + addition)); % opto first  
    vec_SpikesThisTrial_whisker = vec_SpikeTimes(vec_SpikeTimes > vec_StimTimes_stim2(stim) & vec_SpikeTimes < (vec_StimTimes_stim2(stim) + addition));
    vec_SpikesThisTrial_all = vec_SpikeTimes(vec_SpikeTimes > vec_StimTimes_stim1(stim) & vec_SpikeTimes < (vec_StimTimes_stim1(stim) + addition_all)); %look during entire window (both wisker and opto)
    end
    
Spiking_opto(stim) = ~isempty(vec_SpikesThisTrial_opto);
Num_spikes_opto(stim) = numel(vec_SpikesThisTrial_opto);
Spiking_whisker(stim) = ~isempty(vec_SpikesThisTrial_whisker);
Num_spikes_whisker(stim) = numel(vec_SpikesThisTrial_whisker);
Spiking_all(stim) = ~isempty(vec_SpikesThisTrial_all);
Num_spikes_all(stim) = numel(vec_SpikesThisTrial_all);

end

%find the delays between opto and whisker and group them

difference = nan(numStim,1);

for stim2 = 1:numStim
    difference(stim2) = round(vec_StimTimes_stim2(stim2) - vec_StimTimes_stim1(stim2),3);   %delays between stim2 and stim1
end

rounded_difference = round(difference / 0.005) * 0.005;   %round so that it ends in steps of 0.005

[group_idx,Group] = findgroups(rounded_difference);   %split delays in different groups 
numgroups = numel(Group);

% calculate probability of spiking per delay/group

%preallocate
vec_probability_opto = nan(1,numgroups);
vec_probability_whisker = nan(1,numgroups);
vec_probability_all = nan(1,numgroups);
vec_mean_num_spikes_opto = nan(1,numgroups);
vec_mean_num_spikes_whisker = nan(1,numgroups);
vec_mean_num_spikes_all = nan(1,numgroups);
mat_num_spikes = nan(40,numgroups);
mat_contingency_all = nan(6,numgroups);
vec_contingency_opto = nan(5,1);
vec_contingency_whisker = nan(5,1);

for i = 1:numgroups
    groupMembers = find(group_idx == i);
    groupSpikingProbabilities_opto = Spiking_opto(groupMembers);
    groupSpikingProbabilities_whisker = Spiking_whisker(groupMembers);
    groupSpikingProbabilities_all = Spiking_all(groupMembers);
    groupNumSpikes_opto = Num_spikes_opto(groupMembers);
    groupNumSpikes_whisker = Num_spikes_whisker(groupMembers);
    groupNumSpikes_all = Num_spikes_all(groupMembers);
    num_trials = numel(groupNumSpikes_all); 
    % count spikes per trial and accumulate in matrix
    mat_num_spikes(1:num_trials,i) = groupNumSpikes_all';
        % count distribution of spikes number/trial (1 to 5) for statistics
        for spikenum = 0:5
        mat_contingency_all(spikenum+1,i) = numel(find(groupNumSpikes_all==spikenum));
        end

    vec_probability_opto(i) = mean(groupSpikingProbabilities_opto);
    vec_probability_whisker(i) = mean(groupSpikingProbabilities_whisker);
    vec_probability_all(i) = mean(groupSpikingProbabilities_all);
    vec_mean_num_spikes_opto(i) = mean(groupNumSpikes_opto);
    vec_mean_num_spikes_whisker(i) = mean(groupNumSpikes_whisker);
    vec_mean_num_spikes_all(i) = mean(groupNumSpikes_all);
    % do the same for opto/whisker only trials
    if Group(i) == 1
    vec_opto_only_spikes = groupNumSpikes_opto;
    vec_whisker_only_spikes = groupNumSpikes_whisker;
    vec_sum = groupNumSpikes_opto + groupNumSpikes_whisker;
        % also here count distribution of spikes
        for spikenum = 0:5
        vec_contingency_opto(spikenum+1,1) = numel(find(groupNumSpikes_opto==spikenum));  
        vec_contingency_whisker(spikenum+1,1) = numel(find(groupNumSpikes_whisker==spikenum));
        vec_contingency_sum(spikenum+1,1) = numel(find(vec_sum==spikenum));
        end
    end
           
end

%preallocate excel naming
filename = cellstr(files(file).name);
header_one = cellstr('opto_only');
header_two = cellstr('whisker_only');
header_three = cellstr('sum');

%preallocate excel cells to paste contingency data
excel_cell_matrix = strcat('B',num2str(num_trials+5));
excel_cell_opto = strcat('N',num2str(num_trials+5));
excel_cell_whisker = strcat('O',num2str(num_trials+5));
excel_cell_sum = strcat('P',num2str(num_trials+5));

%write excel to current folder
 xlswrite(cellname_excel,filename,file,'A1');
 xlswrite(cellname_excel,Group',file,'B2');
 xlswrite(cellname_excel,mat_num_spikes,file,'B3');
 xlswrite(cellname_excel,header_one,file,'N2');
 xlswrite(cellname_excel,header_two,file,'O2');
 xlswrite(cellname_excel,header_three,file,'P2');
 xlswrite(cellname_excel,vec_opto_only_spikes',file,'N3');
 xlswrite(cellname_excel,vec_whisker_only_spikes',file,'O3');
 xlswrite(cellname_excel,vec_sum',file,'P3');
 xlswrite(cellname_excel,mat_contingency_all,file,excel_cell_matrix);
 xlswrite(cellname_excel,vec_contingency_opto,file,excel_cell_opto);     
 xlswrite(cellname_excel,vec_contingency_whisker,file,excel_cell_whisker);
 xlswrite(cellname_excel,vec_contingency_sum,file,excel_cell_sum);

%make matrices with averaged data 
matrix_probabilities_opto(1:numel(vec_probability_opto),file) = vec_probability_opto';
matrix_probabilities_whisker(1:numel(vec_probability_whisker),file) = vec_probability_whisker';
matrix_probabilities_all(1:numel(vec_probability_all),file) = vec_probability_all';
matrix_spikenumber_opto(1:numgroups,file) = vec_mean_num_spikes_opto';
matrix_spikenumber_whisker(1:numgroups,file) = vec_mean_num_spikes_whisker';
matrix_spikenumber_all(1:numgroups,file) = vec_mean_num_spikes_all';

headers(file) = cellstr(files(file).name);
matrix_groups(1:numel(Group),file) = Group';
end



%% EXPORT TO EXCEL

tdir=uigetdir;  %saving location
cd(tdir)

text1 = cellstr('probabilities opto');
text2 = cellstr('spike number opto');
text3 = cellstr('probabilities whisker');
text4 = cellstr('spike number whisker');
text5 = cellstr('delays');


xlswrite(excelsheet_filename,text1,1,'A2');
xlswrite(excelsheet_filename,text2,1,'A17');
xlswrite(excelsheet_filename,text3,1,'A32');
xlswrite(excelsheet_filename,text4,1,'A46');
xlswrite(excelsheet_filename,text4,1,'A60');
xlswrite(excelsheet_filename,headers,1,'B1');
xlswrite(excelsheet_filename,matrix_probabilities_opto,1,'B2');
xlswrite(excelsheet_filename,matrix_spikenumber_opto,1,'B17');
xlswrite(excelsheet_filename,matrix_probabilities_whisker,1,'B32');
xlswrite(excelsheet_filename,matrix_spikenumber_whisker,1,'B46');
xlswrite(excelsheet_filename,matrix_groups,1,'B60');
xlswrite(excelsheet_filename,headers,2,'B1');
xlswrite(excelsheet_filename,matrix_probabilities_all,2,'B2');
xlswrite(excelsheet_filename,matrix_spikenumber_all,2,'B17');
xlswrite(excelsheet_filename,matrix_groups,2,'B32');



