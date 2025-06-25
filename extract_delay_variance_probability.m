function [Evoked,probability,first_delay,first_variance,mean_number_of_spikes,mean_number_of_spikes_ALLtrials,num_trials,mean_inst_freq,number_of_bursts,ALL_Delays,ALL_first_Delays,Zeta,ZetaP] = extract_delay_variance_probability(vec_SpikeTimes,vec_StimTimes,stimduration,addition,spont_period,intPlot)

% This script extracts delays, variance and probability spiking in response
%to a stimulus
% input: vec_Spiketimes: spike times in s
%        vec_Stimtimes: stim times in s
%        stimduration: duration of stimulus in s
%        addition: window on top of stimulus duration to look for spikes(s)
%        spont_period: window before which spikes are considered
%        spontaneous


 
vec_StimEnd = vec_StimTimes + stimduration;


%Preallocate variables
numstim = numel(vec_StimTimes);
Evoked = cell(numstim,1);  %SpikeTimes of spikes during stim
Delays = cell(numstim,1);  %Delay to onset of spikes during stim 
Spiking = nan(numstim,1);
vec_first_delays = nan(numstim,1);
vec_num_spikes = nan(numstim,1);
vec_inst_freq = nan(numstim,1);

%Check for opto-evoked spikes
for stim = 1:numstim
vec_SpikesThisTrial = vec_SpikeTimes(vec_SpikeTimes > vec_StimTimes(stim) & vec_SpikeTimes < (vec_StimEnd(stim) + addition));
vec_delays = vec_SpikesThisTrial - vec_StimTimes(stim); %delay from stim onset   
Evoked{stim,1} = vec_SpikesThisTrial;
Delays{stim,1} = vec_delays;
Spiking(stim) = ~isempty(vec_SpikesThisTrial);   %spiking or not during stim --> for probability calculation (failures)
vec_num_spikes(stim) = numel(vec_SpikesThisTrial);


    if ~isempty(vec_delays)
    
         %exclude spikes that occur too early --> must be spontaneous
         if vec_delays(1) > spont_period         
         vec_first_delays(stim,1) = vec_delays(1);
         elseif numel(vec_delays) > 1 % then take second spike
         vec_first_delays(stim,1) = vec_delays(2); 
         else                                     %dont count
         vec_first_delays(stim,1) = [nan];   
         end
    
        %calculate instantaneous burst rate if two or more spikes fired
        if numel(vec_delays) > 1
        vec_inst_freq(stim,1) = 1/(vec_delays(2) - vec_delays(1)); %Instantaneous frequency in Hz
        end

    else
    end

end

%Do zeta test

dblBaselineDuration = 0.01;
dblUseMaxDur = 0.05; % full duration (baseline plus during/after stim)
intPlot = 0;

[Zeta,ZetaP,~] = zeta_baseline(dblBaselineDuration,dblUseMaxDur,vec_StimTimes,vec_SpikeTimes,intPlot);

%Cumulate average values for each cell(file) in new vector

probability = numel(find(Spiking ==1))/numstim;
first_delay = (nanmean(vec_first_delays))*1000;  % in ms
first_variance = (nanvar(vec_first_delays*1000)); % in ms2  
mean_number_of_spikes = nanmean(nonzeros(vec_num_spikes)); % if its spiking (not incuding failures!)
mean_number_of_spikes_ALLtrials = nanmean(vec_num_spikes); %including failures
num_trials = numstim;
mean_inst_freq = nanmean(vec_inst_freq);  
number_of_bursts = sum(~isnan(vec_inst_freq));
ALL_Delays = (sort(vertcat(Delays{:})))*1000; % in ms
ALL_Delays (isnan(ALL_Delays)) = []; % remove nans
ALL_first_Delays = (vec_first_delays)*1000; % in ms
ALL_first_Delays(isnan(ALL_first_Delays)) = [];


end