function [sAP,pulse_times_all_recordings] = getpulsetimes(rec_name)
%GETPULSETIMES generates a cell array of pulse times (stumuli) of all recordings
%(opto or whisker (single), or the combination (dual))

% input is rec_name --> looks for a file with that name to open (sAP) in user defined directory, then detects
% all opto recordings ('RunOptoStim'), loops through all recordings and
% generates a list of all stimulus ON times based on the info of the
% stimuli (at which frequencies + durations pulses were delivered), saves
% them in a cell array "pulse_times_all_recordings" (output of function)

file_path = strcat('\\vs03.herseninstituut.knaw.nl\VS03-AXS-1\Jamann_Longrange_projections_1\Analysis\Neuropixels\RecordingProcessor\new synthesis with bombcell\Ctrl');
load(rec_name);

%%
number_of_recordings = numel(sAP.cellBlock); %all recordings (including other types of stimuli)
 
for recidx = 1:number_of_recordings
rec_types(recidx) = cellstr(sAP.cellBlock{1, recidx}.strExpType);   %looking for opto stim or optowhisker stim
single_rec_logical(recidx) = contains(rec_types(recidx),'RunOptoStim');
dual_rec_logical(recidx) = contains(rec_types(recidx),'RunOptoWhiskerStim');
end

single_recs = find(single_rec_logical(:) == 1);  % list of all opto/whisker only recordings
number_of_single_recordings = numel(single_recs); %number of single opto/whisker recordings
dual_recs = find(dual_rec_logical(:) == 1);  % list of all dual (opto + whisker) recordings
number_of_dual_recordings = numel(dual_recs);  %number of dual recordings

pulse_times_all_recordings = cell(1,number_of_recordings);
sAP = sAP;

%% Opto/whisker recordings
%loop runs through each opto/whisker only recording and concatenates pulse times
%according to the info on stimuli in cellBlock and adds those times to the
%NI timestamps ("vecStimonTime")

for singleInt = 1:number_of_single_recordings

% accumulate all info for this recording from sAP.cellBlock
single_rec_idx = single_recs(singleInt);
number_of_trials = numel(sAP.sSources.cellBlock{1, single_rec_idx}.structEP.vecStimOnTime);
pulse_dur_single = sAP.cellBlock{1, single_rec_idx}.vecPulseDur(1);  %pulse duration (in s) only if same pulse length thoughout!
ITI_all_trials = sAP.cellBlock{1, single_rec_idx}.cellPulseITI; %interval (between pulses) for all trials
num_pulse = sAP.cellBlock{1, single_rec_idx}.intRepsPerPulse;   % number of pulses in one pulse train
num_stim = sAP.cellBlock{1, single_rec_idx}.intStimTypes;    %number of stimulation frequencies

% CONCATENATE ALL OPTO STIMULATION TIMES FROM THIS RECORDING
%find all stimulus ON times for one trial from the pulse infos

abs_all_pulse_times = zeros(num_stim*num_pulse, number_of_trials); %preallocate matrix for all trials

for trial_idx = 1:number_of_trials
    
    ITI_info = ITI_all_trials{trial_idx};  %order of frequencies within that trial (are shuffled randomly)
    all_trains = zeros(num_pulse,num_stim);  %matrix of pulse onsets for all trains within one trial (seconds)
    vecTrain = zeros(1,num_pulse); %one pulse train (e.g. five opto pulses per train)
    
    %for each pulse onset we need to add the info on the pulse frequency,
    %pulse duration to the previous pulse train: always 1 s in between trials
    %plus one ITI and one pulse duration
    
    for freq_idx = 1:num_stim    %freq_idx = the n different frequencies
        
        if freq_idx == 1
             for a = 2:num_pulse
                vecTrain(a) = (a-1)*ITI_info(freq_idx) + (pulse_dur_single)*(a-1);
                all_trains(:,freq_idx) = vecTrain;   
            end
        else
            
            for a = 2:num_pulse
                onset = all_trains(num_pulse,freq_idx-1)+ sAP.cellBlock{1, single_rec_idx}.dblPulseWaitSignal + pulse_dur_single + ITI_info(freq_idx-1);
                vecTrain(a) = onset + (a-1)*ITI_info(freq_idx) + (pulse_dur_single)*(a-1);
                vecTrain(1) = onset;
                all_trains(:,freq_idx) = vecTrain;
            
            end
       
        end
    end
   
%here we add the pulsetimes to the trial-onset according to the NI timestamp
%note: if trials were aborted earlier, then the last trial onset time is not saved, but we can 
%deduct it from the difference in onset time between trial
            
            format long  %this is to prevent rounding
            %vecStimonTime = (sAP.cellBlock{1, single_rec_idx}.ActOnNI(trial_idx)-sAP.cellBlock{1, 1}.T0);
            vecStimonTime = (sAP.cellBlock{1, single_rec_idx}.vecStimOnTime(trial_idx));
            trial_onset = vecStimonTime; %true trial onset 
            all_pulse_times = all_trains(:);
            format long
            abs_all_pulse_times(:,trial_idx) = all_pulse_times + trial_onset; 
            
            

end
    
    all_pulsetimes_combined_rec = (abs_all_pulse_times(:));   %make list of all pulse times in one rec
    pulse_times_all_recordings{1,single_rec_idx} = all_pulsetimes_combined_rec;  %store in cell array --> output of function
    
end
    
%% Dual recordings
%loop runs through each dual recording and concatenates the pulse times
%(first stimulus) according to the info on stimuli in cellBlock and adds those times to the
%NI timestamps ("vecStimonTime")

for dualInt = 1:number_of_dual_recordings

% accumulate all info for this recording from sAP.cellBlock
dual_rec_idx = dual_recs(dualInt);
number_of_dual_trials = numel(sAP.sSources.cellBlock{1, dual_rec_idx}.structEP.vecStimOnTime);
pulse_dur_dual = sAP.cellBlock{1, dual_rec_idx}.vecPulseDur_opto(1);  %pulse duration (in s) 
ITI_dual = sAP.cellBlock{1, dual_rec_idx}.vecPulseITI + sAP.cellBlock{1, dual_rec_idx}.dblPulseWaitSignal + pulse_dur_dual; %interval between pulses (ITI, pulsewait and pulseduration)
num_delays = numel(sAP.cellBlock{1, dual_rec_idx}.cellPulseDelay{1, 1});

% CONCATENATE ALL STIMULATION TIMES FROM THIS RECORDING
%preallocate matrix for all trials and delays 
all_trains_dual = nan(num_delays,number_of_dual_trials);
stim_on = nan(1,num_delays);

    for delay_idx = 1:num_delays
    stim_on(delay_idx) = (delay_idx-1) * ITI_dual;
    end
 
for trial_idx = 1:number_of_dual_trials
    
        trial_onset = (sAP.cellBlock{1, dual_rec_idx}.vecStimOnTime(trial_idx)); %true trial onset        
        all_trains_dual(:,trial_idx) = trial_onset + stim_on; 
end
        all_pulsetimes_combined_dual = (all_trains_dual(:));
        pulse_times_all_recordings{1,dual_rec_idx} = all_pulsetimes_combined_dual;
end
end

