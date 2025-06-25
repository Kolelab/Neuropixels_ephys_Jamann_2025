%%LOAD DATA

%User locates target directory with exel files
tdir=uigetdir;
cd(tdir)

files = dir('*.xlsx');

filenumber = numel(files);
stimduration = 0.040; %indicate length of opto stim in s plus 5ms or so since sometimes delay is longert
% min_delay = 0.05; %to exclude random spikes that appear too early

%% EDIT TABLE
%load datatable, reassign headers, make new columns with pulse and
%spiketimes (because at each episode latency starts at 0 again)

oneHz = zeros(5,filenumber);
tenHz = zeros(5,filenumber);
twentyHz = zeros(5,filenumber);
fiftyHz = zeros(5,filenumber);
hundredHz = zeros(5,filenumber);
mean_delays = zeros(1,filenumber);
variance_delays = zeros(1,filenumber);
probability_first_pulse = zeros(1,filenumber);


for file = 1:filenumber
     
T = readtable(files(file).name);     %generatetable with data of file a
T.Properties.VariableNames{4} = 'PulseTimes';

T.PulseTimes = (T.Latency_s_ + (T.Episode_-1)*20);   %concatenate all pulse and spike times relative to 0
T.SpikeTimes = (T.Latency_s__1 + (T.Episode__1-1)*20);
T.PulseEnd = T.PulseTimes + stimduration;

SpikeTimes = T.SpikeTimes;
%% Check for opto-evoked spikes

T.Spiking(:) = 0;

for i = 1:numel(T.Spiking)
    
T.Spiking(i) = any(SpikeTimes > T.PulseTimes(i) & SpikeTimes < T.PulseEnd(i));  %finds any pulse where a spike followed 

end

%% IDENTIFY FREQUENCIES OF STIMULATION
%based on the interspike interval (ISI) between pulses
%Check table afterwards for mistakes in event detection!

for n = 1:(numel(T.PulseTimes)-1)
    
T.ISI(n) = round(T.PulseTimes(n+1)-T.PulseTimes(n)-0.02,4);


if T.ISI(n) <= 1
    
T.Frequency(n) = round(1/(T.ISI(n)));

else
    if n == 1
        T.Frequency(n) = 0;
    else
    T.Frequency(n) = T.Frequency(n-1);
    end
end
    

end

%% SPIKING PROBABILITY
%%loops through all frequencies, finds for pulse 1 to 5 if spike was fired
%%or not and calculates overall probability

probabilities = zeros(5);
frequencies = [1,10,20,50,100];  %edit frequencies if different ones were used

for o = frequencies

cell.stim = find(T.Frequency == o);  %create structure cell, find all stimuli of that frequency
cell.firing = T.Spiking(cell.stim);  %for that pulse, was it firing or not
cell.pulse = zeros(numel(cell.stim),1);  


for m = 1:(numel(cell.stim)-5)     %loop to find pulse number 1 to 5 --> if 5 consecutive times the same index
    
    if m ==1 
        cell.pulse(m) = 1;
    else
    if (cell.stim(m+4)-cell.stim(m)== 4) && (cell.stim(m+5) - cell.stim(m)> 5);
        cell.pulse(m)= 1;
    else cell.pulse(m) = cell.pulse(m-1)+1;
    end
    end
end

%%find the number of spikes evoked by each pulse

for k = 1:5

pulse_number = find(cell.pulse == k);  %index of pulses 1 to 5 of that specific frequency
spikes = find(cell.firing(pulse_number) == 1); %how many of the times was a spike fired  
probability(k) = numel(spikes)/numel(pulse_number); %probability is number of spikes divided by number of total pulses

end

idx_frequency = find(frequencies == o);
probabilities(idx_frequency,:) = probability;

end

%creates matrixes for each frequency and fill in data from each cell

oneHz(:,file) = probabilities(1,:); 
tenHz(:,file) = probabilities(2,:);
twentyHz(:,file) = probabilities(3,:);
fiftyHz(:,file) = probabilities(4,:);
hundredHz(:,file) = probabilities(5,:);
headers(file) = cellstr(files(file).name);

%%find probability of all first pulses

for z = 2:numel(T.Frequency)
    
if T.Frequency(z) == T.Frequency(z-1)
    T.FirstPulse(z) = 0;
else
    T.FirstPulse(z) = 1;
    
end
end

first_pulses = find(T.FirstPulse == 1);
spiking_first_pulse = find(T.Spiking(first_pulses)==1);
probability_first_pulse(file) = numel(spiking_first_pulse)/(numel(first_pulses));


%% DELAYS & VARIANCE
%find every spike during pulse and calculate delay, from that calculate mean delay & variance

optopulses = find(T.Spiking == 1);%index of all pulses that evoked a spike
optopulsetimes = T.PulseTimes(optopulses);
spiketimes_list = T.SpikeTimes;
delays = zeros(numel(optopulses),1);

for p = 1:numel(optopulses)
    
evoked_spike_index = find(spiketimes_list > optopulsetimes(p),1);  %find the first spike after each pulse onset
delays(p) = spiketimes_list(evoked_spike_index) - optopulsetimes(p);
mean_delays(file) = mean(delays)*1000;  %calculates mean of all delays of one neuron
variance_delays(file) = var(delays)*1000; %calculates variance of all delays of one neuron
delaycells.all{file} = (delays)*1000;
relative_delays = delays/(mean_delays(file))*1000;  %each delay relative to mean --> to be able to plot histogram over multiple cells
delaycells.relative{file} = relative_delays;


end


 
end



%% SAVE EXCEL SHEET
%Save output to excel sheet (created in the beginning)

tdir=uigetdir;  %saving location
cd(tdir)
pulses = (1:5)';
text1 = cellstr('average delay');
text2 = cellstr('variance');
text3 = cellstr('all delays');
text4 = cellstr('relative delays');

excelsheet_filename = 'firing_probability4.xlsx';
xlswrite(excelsheet_filename,headers,1,'B1');
xlswrite(excelsheet_filename,oneHz,1,'B3');
xlswrite(excelsheet_filename,tenHz,1,'B10');
xlswrite(excelsheet_filename,twentyHz,1,'B17');
xlswrite(excelsheet_filename,fiftyHz,1,'B24');
xlswrite(excelsheet_filename,hundredHz,1,'B31');
xlswrite(excelsheet_filename,pulses,1,'A3');
xlswrite(excelsheet_filename,pulses,1,'A10');
xlswrite(excelsheet_filename,pulses,1,'A17');
xlswrite(excelsheet_filename,pulses,1,'A24');
xlswrite(excelsheet_filename,pulses,1,'A31');
xlswrite(excelsheet_filename,probability_first_pulse,1,'B38');
xlswrite(excelsheet_filename,text1,2,'A2');
xlswrite(excelsheet_filename,text2,2,'A3');
xlswrite(excelsheet_filename,text3,2,'A5');
xlswrite(excelsheet_filename,text4,3,'A5');
xlswrite(excelsheet_filename,headers,2,'B1');
xlswrite(excelsheet_filename,mean_delays,2,'B2');
xlswrite(excelsheet_filename,variance_delays,2,'B3');
xlswrite(excelsheet_filename,headers,3,'B1');


% letters_vec = ['B5'; 'C5'; 'D5'; 'E5'; 'F5'; 'G5'; 'H5'; 'I5'; 'J5'; 'K5'; 'L5'; 'M5'; 'N5'; 'O5'; 'P5'; 'Q5'; 'R5', 'S5', 'T5', 'U5', 'V5', 'W5', 'X5', 'Y5']
% letters_str = cellstr(letters_vec);
% 
% for q = 1:numel(files)
%     
% xlswrite(excelsheet_filename,delaycells.all{q},2,letters_str{q});
% xlswrite(excelsheet_filename,delaycells.relative{q},3,letters_str{q});
% 
% end




