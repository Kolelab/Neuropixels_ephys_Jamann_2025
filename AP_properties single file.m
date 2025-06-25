%% Select file
clear
dir = uigetdir;
cd(dir)

axo = load('210225_cell1 003 Copy Export.mat');                %add file name

clear('dir')



% set parameters                   

Fs = inputdlg('Sample Rate (in Hz):');
Fs = str2double(Fs);
repetitions = inputdlg('Number of Repetitions:');
repetitions = str2double(repetitions);
CI_onset = inputdlg('current injection onset (s):');
CI_onset = str2double(CI_onset);
CI_onset = CI_onset*Fs;

fns = fieldnames(axo);
i=0;
index = 2;                          %the initial index refers to the first file in the "axo" struct that contains meaningful data !CHANGE IF NEEDED!
while i < repetitions               %loop to obtain data from axo struct
    i=i+1;
    axo_data{i} = axo.(fns{index});
    index = index + 2;              %index is by default set to increase by 2, as two graphs are generated per sweep of which only one is useful
end
clear('index','i')


% AP properties                     %ignore Warning messages

%peak properties
j=0;
all_peaks = zeros(1,repetitions);
all_locations = zeros(1,repetitions);
while j < repetitions
    j = j+1;
    [pks,locs] = findpeaks(axo_data{1,j},'MinPeakHeight',0);
    if pks > 0
    else
        pks = 0;
        locs = 0;
    end
    all_peaks(j) = pks;
    all_locations(j) = locs/Fs;
end
clear('j','pks','locs')
peaks = nonzeros(all_peaks);
locations = nonzeros(all_locations);


%rise time

k = 0;
while k < length(locations)
    k = k+1;
    deriv = diff(axo_data{k})*Fs;
    seekloc = locations(k)*Fs-(0.002*Fs);
    end_rise = find(deriv(seekloc:locations(k)*Fs)>50,1);
    AP_start = (locations(k)*Fs - 0.002*Fs) + end_rise;
    rise_time = AP_start - CI_onset;
    rise_time = rise_time/Fs;
    rise_times(k) = rise_time;
end
clear('k','AP_start','rise_time','end_rise','seekloc','deriv')        %remove temp. variables


%amplitude and half-width

l = 0;
m = 0;
while l < length(axo_data)
    l = l+1;
    [pks,locs] = findpeaks(axo_data{1,l},'MinPeakHeight',0);
    if pks > 0
    else
        pks = 0;
        locs = 0;
    end
    if pks == 0
    else
        
        %amplitude
        m = m+1;
        AP_start = CI_onset + rise_times(m)*Fs;
        AP_treshold = axo_data{l}(AP_start);
        amplitude(m) = peaks(m) - AP_treshold;
        
        
        %half-width
        half_amp = AP_treshold + (amplitude(m)/2);
        half_ap_rise = find(axo_data{l} >= half_amp,1);
        half_ap_down = find(axo_data{l}(locations(m)*Fs:length(axo_data{l}))>= half_amp,1) + (locations(m)*Fs);
        half_width(m) = (half_ap_down - half_ap_rise)/Fs;
        
    end
end
clear('m','l','locs','pks')

% data summary

success_rate = length(locations) / length(all_locations);
mean_amplitude = mean(amplitude);
mean_half_width = mean(half_width);
mean_rise_time = mean(rise_times);


data_summary = struct("AP_success_rate", success_rate, "mean_amplitude", mean_amplitude, ...
                    "mean_halfwidth", mean_half_width, "mean_rise_time", mean_rise_time, ...
                    "all_amplitudes", amplitude, "all_halfwidths", half_width, ...
                    "all_rise_times_", rise_times, "AP_peak_times", locations);
display(data_summary)     







