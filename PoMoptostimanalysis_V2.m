
% This script analyses exported data from axograph saved as .xlsx files.
% One column has the detected pulses,one the detected spikes
% It will output average delays, variances and probabilities of spiking and burst per cell
% The results can be saved in a summary excel file, separated by file


%% LOAD DATA

% User locates target directory with exel files
tdir = uigetdir;

cd(tdir)

files = dir('*.xlsx');

filenumber = numel(files);
stimduration = 0.020; %indicate length of opto stim (in s)
excelsheet_filename = 'delays_variance_probality_PoM_Ctrl_240617';
excelsheet_filename_two = 'bursts_Ctrl_added_240617_3';
%offset_window = stimduration + addition;
prompt = 'Analyse bursts? Yes = 1, No = 0';
analyse_burst = cell2mat(inputdlg(prompt));


%% Load spike/stimtimes

[SpikeTimes,StimTimes,totaltime] = find_optostim(files,stimduration);

%% Extract data
    
for file = filenumber

[basePath, filename, ~] = fileparts(files(file).name); 
%extract spike and stimtimes
vec_SpikeTimes = cell2mat(SpikeTimes(1,file));  
vec_StimTimes = cell2mat(StimTimes(1,file));  

dblBaselineDuration = 0.02;
% plotRasterdots(vec_SpikeTimes,vec_StimTimes,0.05,dblBaselineDuration);
% filename_raster = strcat(filename,' raster');
% current_xticks = get(gca, 'XTick');
% new_xticks = current_xticks - 0.02;
% %set(gca, 'XTick', new_xticks);
% new_xtick_labels = arrayfun(@(v) sprintf('%.2f', v), new_xticks, 'UniformOutput', false);
% set(gca, 'XTickLabel', new_xtick_labels)
% savefig(filename_raster)
% close

% psth = mpsth(vec_SpikeTimes,vec_StimTimes,'fr',0,'binsz',1,'pre',20,'post',30,'chart',1);
% filename_psth = strcat(filename,' psth');
% savefig(filename_psth)
% close


stimduration = 0.020; %indicate length of opto stim (in s)
addition = 0.02; %indicate the whole window to look for spikes (on top of the stim length) (in s)
end_window = stimduration + addition;
spont_period = 0.003; % if earlier than this period after onset then spontaneous


%function to extract delay, variance, probability of spiking
[Evoked,probability,first_delay,first_variance,mean_number_of_spikes,mean_number_of_spikes_ALLtrials,num_trials,mean_inst_freq,number_of_bursts,ALL_Delays,ALL_first_Delays,Zeta,ZetaP] = extract_delay_variance_probability(vec_SpikeTimes,vec_StimTimes,stimduration,addition,spont_period);


%Cumulate average values for each cell (file) into new vector
vec_probability(file) = probability;
vec_mean_delay(file) = first_delay;
vec_variance_delay(file)= first_variance;
vec_mean_number_of_spikes(file) = mean_number_of_spikes;
vec_mean_number_of_spikes_ALLtrials(file) = mean_number_of_spikes_ALLtrials;
vec_num_trials(file) = num_trials;
headers(file) = cellstr(filename);
vec_mean_inst_freq(file) = mean_inst_freq;  
vec_number_of_bursts(file) = number_of_bursts;
vec_zeta(file) = Zeta;
vec_zetaP(file) = ZetaP;
ALL_Delays_cat(1:numel(ALL_Delays),file) = ALL_Delays';
ALL_first_Delays_cat(1:numel(ALL_first_Delays),file) = ALL_first_Delays';

if analyse_burst == '1'

%function to extract bursting properties from spike and stimtimes
[start_first,end_first,start_second,end_second,mean_delay_first,mean_variance_first,mean_delay_second,mean_variance_second,probability_first,probability_second,burst_number,burst_probability,frequency_first_two_spikes] = extract_bursts(filename,vec_SpikeTimes,vec_StimTimes);

%Cumulate average burst values
vec_delay_first(file) = mean_delay_first;
vec_variance_first(file) = mean_variance_first;
vec_delay_second(file) = mean_delay_second;
vec_variance_second(file) = mean_variance_second;
vec_probability_first(file) = probability_first;
vec_probability_second(file) = probability_second;
vec_probability_burst(file) = burst_probability;
vec_burst_number(file) = burst_number; %this is only the first and second spike in window, above its all, so also later spikes!
vec_burst_frequency_first_two(file) = frequency_first_two_spikes;
vec_start_first(file) = start_first;
vec_end_first(file) = end_first;
vec_start_second(file) = start_second;
vec_end_second(file) = end_second;

else 
end

% spontaneous bursts

Evoked_spikes = sort(vertcat(Evoked{:}));
vec_spontaneous = setdiff(vec_SpikeTimes,Evoked_spikes);

num_spont = numel(vec_spontaneous);

if ~isempty(vec_spontaneous)
    
% Parameters
spike_threshold = 50; % Minimum spike frequency to be considered a burst (in Hz)
%refractory_period = 0.01; % Minimum time gap between bursts (in seconds)

% Calculate spike intervals
spike_intervals = diff(vec_spontaneous);
spike_frequencies = 1 ./ spike_intervals;

% Find the indices where the spike frequency exceeds the threshold
burst_indices = find(spike_frequencies > spike_threshold);
differences = diff(burst_indices);
burst_remove_indices = (burst_indices(find(differences == 1)))+1; %remove values that belong to a burst with several APs
burst_indices_removed = burst_indices(~ismember(burst_indices,burst_remove_indices));
Instfreq_spont_bursts = spike_frequencies(burst_indices_removed); 
vec_Instfreq_spont_bursts(file) = mean(Instfreq_spont_bursts);
num_spont_bursts(file) = numel(Instfreq_spont_bursts); 
vec_burstpermin(file) =  numel(Instfreq_spont_bursts)/totaltime;

else
end

end






%% EXCTRACT all dalays of all spikes
% 
vec_ALL_Delays = sort(ALL_Delays_cat(:));
vec_ALL_Delays = vec_ALL_Delays(vec_ALL_Delays > 0);
vec_ALL_first_delays = sort(ALL_first_Delays_cat(:));
vec_ALL_first_delays = vec_ALL_first_delays(vec_ALL_first_delays>0);


%% SAVE EXCEL SHEET ONE
%Save output to excel sheet (created in the beginning)

tdir=uigetdir;  %saving location
cd(tdir)
text1 = cellstr('mean delay');
text2 = cellstr('variance');
text3 = cellstr('mean number of spikes');
text4 = cellstr('mean number of spikes all trials');
text5 = cellstr('probability');
text6 = cellstr('number of trials');
text7 = cellstr('burst frequency');
text8 = cellstr('burst number');
text9 = cellstr('spont burst number');
text10 = cellstr('spont burst freqeuncy');
text11 = cellstr('spont burst per min');
text12 = cellstr('zeta');
text13 = cellstr('zeta P');
text14 = cellstr('start window');
text15 = cellstr('end window');


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
xlswrite(excelsheet_filename,text12,1,'A13');
xlswrite(excelsheet_filename,text13,1,'A14');
xlswrite(excelsheet_filename,text14,1,'A15');
xlswrite(excelsheet_filename,text14,1,'A16');

xlswrite(excelsheet_filename,headers,1,'B1');
xlswrite(excelsheet_filename,vec_mean_delay,1,'B2');
xlswrite(excelsheet_filename,vec_variance_delay,1,'B3');
xlswrite(excelsheet_filename,vec_mean_number_of_spikes,1,'B4');
xlswrite(excelsheet_filename,vec_mean_number_of_spikes_ALLtrials,1,'B5');
xlswrite(excelsheet_filename,vec_probability,1,'B6');
xlswrite(excelsheet_filename,vec_num_trials,1,'B7');
xlswrite(excelsheet_filename,vec_mean_inst_freq,1,'B8');
xlswrite(excelsheet_filename,vec_number_of_bursts,1,'B9');
xlswrite(excelsheet_filename,num_spont_bursts,1,'B10');
xlswrite(excelsheet_filename,vec_Instfreq_spont_bursts,1,'B11');
xlswrite(excelsheet_filename,vec_burstpermin,1,'B12');
xlswrite(excelsheet_filename,vec_zeta,1,'B13');
xlswrite(excelsheet_filename,vec_zetaP,1,'B14');
xlswrite(excelsheet_filename,spont_period,1,'B15');
xlswrite(excelsheet_filename,end_window,1,'B16');

%% SAVE excel sheet 2


tdir=uigetdir;  %saving location
cd(tdir)

text1 = cellstr('delay first spike');
text2 = cellstr('variance first spike');
text3 = cellstr('probability first spike');
text4 = cellstr('delay second spike');
text5 = cellstr('variance second spike');
text6 = cellstr('probability second spike');
text7 = cellstr('probability bursting');
text8 = cellstr('number of bursts');
text9 = cellstr('frequency between first and second');
text10 = cellstr('window first');
text11 = cellstr('window second');

xlswrite(excelsheet_filename_two,text1,1,'A2');
xlswrite(excelsheet_filename_two,text2,1,'A3');
xlswrite(excelsheet_filename_two,text3,1,'A4');
xlswrite(excelsheet_filename_two,text4,1,'A5');
xlswrite(excelsheet_filename_two,text5,1,'A6');
xlswrite(excelsheet_filename_two,text6,1,'A7');
xlswrite(excelsheet_filename_two,text7,1,'A8');
xlswrite(excelsheet_filename_two,text8,1,'A9');
xlswrite(excelsheet_filename_two,text9,1,'A10');
xlswrite(excelsheet_filename_two,text10,1,'A11');
xlswrite(excelsheet_filename_two,text11,1,'A13');
xlswrite(excelsheet_filename_two,headers,1,'B1');
xlswrite(excelsheet_filename_two,vec_delay_first,1,'B2');
xlswrite(excelsheet_filename_two,vec_variance_first,1,'B3');
xlswrite(excelsheet_filename_two,vec_probability_first,1,'B4');
xlswrite(excelsheet_filename_two,vec_delay_second,1,'B5');
xlswrite(excelsheet_filename_two,vec_variance_second,1,'B6');
xlswrite(excelsheet_filename_two,vec_probability_second,1,'B7');
xlswrite(excelsheet_filename_two,vec_probability_burst,1,'B8');
xlswrite(excelsheet_filename_two,vec_burst_number,1,'B9');
xlswrite(excelsheet_filename_two,vec_burst_frequency_first_two,1,'B10');
xlswrite(excelsheet_filename_two,vec_start_first,1,'B11');
xlswrite(excelsheet_filename_two,vec_end_first,1,'B12');
xlswrite(excelsheet_filename_two,vec_start_second,1,'B13');
xlswrite(excelsheet_filename_two,vec_end_second,1,'B14');



