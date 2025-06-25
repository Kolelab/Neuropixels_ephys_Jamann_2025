function [start_first,end_first,start_second,end_second,mean_delay_first,mean_variance_first,mean_delay_second,mean_variance_second,probability_first,probability_second,burst_number,burst_probability,frequency_first_two_spikes] = extract_bursts(filename,vec_SpikeTimes,vec_StimTimes)

% This script extracts the bursts in response to a stimulus
%to a stimulus
% input: vec_Spiketimes: spike times in s
%        vec_Stimtimes: stim times in s
%user input: wndows to look for first and second spike based on raw data
% output: delay, variance, probability of first and second spikes



%ask user for input on search windows for first and second spike
answer1 = inputdlg(strcat("start window first spike ",filename));
answer2 = inputdlg(strcat("end window first spike ",filename));
answer3 = inputdlg(strcat("start window second spike ",filename));
answer4 = inputdlg(strcat("end window second spike ",filename));
start_first = (str2double(answer1{1}))/1000; %in seconds!
end_first = (str2double(answer2{1}))/1000;
start_second = (str2double(answer3{1}))/1000;
end_second = (str2double(answer4{1}))/1000;


%Preallocate variables
numstim = numel(vec_StimTimes);
vec_first_delays_burst = nan(numstim,1);
vec_second_delays_burst = nan(numstim,1);
vec_bursting = zeros(numstim,1);
vec_frequency = nan(numstim,1);

%Check for opto-evoked bursts: first or second spike transmitted?
for stim = 1:numstim
    start_stim = vec_StimTimes(stim);
    FirstSpike_thistrial = vec_SpikeTimes(vec_SpikeTimes < (start_stim +end_first) & vec_SpikeTimes > (start_stim +start_first));  %should be always just one!
    SecondSpike_thistrial = vec_SpikeTimes(vec_SpikeTimes < (start_stim + end_second) & vec_SpikeTimes > (start_stim + start_second));

    if ~isempty(FirstSpike_thistrial)
    vec_first_delays_burst(stim) = (FirstSpike_thistrial - start_stim)*1000; %delay from stim onset in ms 
    else
    end
  
    if ~isempty(SecondSpike_thistrial)
    vec_second_delays_burst(stim) = (SecondSpike_thistrial - start_stim)*1000; %delay from stim onset in ms 
    else
    end

    if ~isempty(FirstSpike_thistrial) && ~isempty(SecondSpike_thistrial) %if two spikes fired
        delay = SecondSpike_thistrial*1000 - FirstSpike_thistrial*1000; 
        vec_frequency(stim) = (1/delay)*1000; %inst frequency between one and two
        if delay < 10 %only if above 100 Hz is burst!
        vec_bursting(stim) = 1;   
        else
        end
    end

end

mean_delay_first = nanmean(vec_first_delays_burst);
mean_variance_first = nanvar(vec_first_delays_burst);
mean_delay_second = nanmean(vec_second_delays_burst);
mean_variance_second = nanvar(vec_second_delays_burst);
probability_first = ((sum(~isnan(vec_first_delays_burst)))/numstim)*100;
probability_second = ((sum(~isnan(vec_second_delays_burst)))/numstim)*100;
burst_number = sum(vec_bursting);
burst_probability = (burst_number/numstim)*100;
frequency_first_two_spikes = nanmean(vec_frequency);


end