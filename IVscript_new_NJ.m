%% IMPORT DATA

%Locate target directory

tdir=uigetdir;
cd(tdir)

% Organize files
fname = dir('*.axtx');

% Convert all .axtx files to .txt files so that they can be read by Matlab
for file_c = 1:numel(fname) 
    
   info = fname(file_c).name;
   text = strrep(info,'.axtx','.txt');
   copyfile(info,text); 
end


%% READ DATA AND EXTRACT PARAMETERS

fname_txt = dir('*.txt');
filenum = numel(fname_txt);
stim_dur = 0.5;   % duration of current injection, adjust if longer!
current_steps = [-200:50:1250];    %modify if different step size/number

%preallocate matrices for parameters

 numAPmatrix = zeros(60,filenum);                     % Total numer of APs
 Frequencies = zeros(60,filenum);                % firing frequency
 ISIavg = zeros(60,filenum);                     % Average of all ISIs in a given Episode
 ISImin = zeros(60,filenum);                     % Maximum value of all ISIs in a given Episode
 MaxIF = zeros(60,filenum);                      % Maximum Instantaneous Frequency in a given Episode
 MeanFrequency = zeros(60,filenum);              % mean frequency based on average ISI
 vec_rheobase = zeros(1,filenum);                  %rheboase current
 vec_rheobase_frequency = zeros(1,filenum);         %average frequency at rheobase current
 mat_fmax = zeros(30,filenum);
 mat_fmax_avg = zeros(30,filenum);
 
for file = 1:filenum 
   
   spiking_eps = [];
   filename = fname_txt(file).name;
   fid = fopen(filename,'r');                       % read the .txt file
   Cell_data = textscan(fid,'%s %s %s','delimiter','\t') ;  % load text file into cell array
   fclose(fid) ;                                    % close the .txt file
   idx_episode = find(contains(Cell_data{1},'Episode'));              % find the rows that have 'Episode' in them    
   lastep = idx_episode(end)+2;
   epnames = Cell_data{1}(idx_episode);                                % Names/Numbers of the Episodes
   epnumber = numel(epnames);
   vec_spiking = zeros(1,epnumber);
   
   
 % Calculating parameters from list
 
 for episode = 1:epnumber                    
            vec_ISI = [];
            idx_startep = idx_episode(episode) +2;                %actual values start two cells down
            if episode < epnumber
            idx_endep = idx_episode(episode+1) -1;
            else 
            endep = numel(Cell_data{1,2});
            end
            vec_spiketimes = str2double(Cell_data{1,2}(idx_startep:idx_endep));    % Select the values for this episode only
            numAPs = numel(vec_spiketimes);
            
            if numAPs == 0                            % Ask whether timepoints is empty, which would mean there were no action potentials
            continue                                  % If timepoints is not empty, it contains an actual value; If timepoints is empty, skip to the next iteration
            
             elseif numAPs == 1          % If there is only 1 AP, there can'b be any ISI/frequency calculations
             vec_spiking(1,episode) = 1;    %logical array if spiking or not --> to find rheobase
             numAPmatrix(episode,file) = numAPs;
             Frequencies(episode,file) = numAPs/stim_dur;
             ISIavg(episode,file) = 0;                                                   
             ISImin(episode,file) = 0;                                                          
             MaxIF(episode,file) = 0;  
             MeanFrequency(episode,file) = 2;
                                        
            else     
             
             vec_spiking(1,episode) = 1; %logical array if spiking or not --> to find rheobase
             for idx_AP = 2:numAPs                        % Calculate the ISIs
             vec_ISI(idx_AP-1) = vec_spiketimes(idx_AP) - vec_spiketimes((idx_AP-1));   % Calculate ISIs
             end
             
             numAPmatrix(episode,file) = numAPs;
             Frequencies(episode,file) = numAPs/stim_dur;
             ISIavg(episode,file) = mean(vec_ISI);                                                   
             ISImin(episode,file) = min(vec_ISI);                                                         
             MaxIF(episode,file) = 1/max(vec_ISI);  
             MeanFrequency(episode,file) = 1/mean(vec_ISI);   
            
            end
    
             
 end

 %find rheobase
            spiking_eps = find(vec_spiking == 1);  % find the spiking episodes
            vec_rheobase(1,file) = current_steps(spiking_eps(1));  %find the first spiking episode
            vec_rheobase_frequency(1,file) = MeanFrequency(spiking_eps(1),file); %frequency at rheobase
            
%calculate max f' (Hz/pA)
            for ep = 2:epnumber
            mat_fmax(ep,file) = (Frequencies(ep,file)- Frequencies(ep-1,file))/50;   %50pA current steps
            mat_fmax_avg(ep,file) = (MeanFrequency(ep,file)- MeanFrequency(ep-1,file))/50;
            end 
            [fmax,fmax_idx] = max(mat_fmax(1:15,file));
            vec_fmax(1,file) = fmax(1);                             %max increase Hz/pA (fmax)
            vec_current_fmax(1,file) = current_steps(fmax_idx(1));  %current at the max increase            
 
end



%% Create EXCEL

%write excel sheet with seperate sheets for each parameter, the headers are
%the filenames, rows are input steps in pA, top left shows which parameter

tdir=uigetdir;
cd(tdir)

excelsheet_filename = 'all_data_DREADD_230413.xlsx';

headers = {};
for header_idx = 1:filenum
    headers{header_idx} = fname_txt(header_idx).name;
end

rownames = [-200:50:1250 -200:50:1250]';%edit if you have different pA steps
parameter = ['APnumber','Frequency','Average ISI','Min ISI','Max ISI','Max Insst Frequency'];

for  g = 1:numel(parameter)
    
    writematrix(parameter(g),excelsheet_filename,'Sheet',g,'Range','A1');
    writecell(headers,excelsheet_filename,'Sheet',g,'Range','B1:GG1');
    writematrix(rownames,excelsheet_filename,'Sheet',g,'Range','A2:A62');%edit if you have more steps/reps
    
end

writematrix(numAPmatrix,excelsheet_filename,'Sheet',1,'Range','B2:GG62');%edit if you have more steps/reps
writematrix(Frequencies,excelsheet_filename,'Sheet',2,'Range','B2:GG62');%edit if you have more steps/reps
writematrix(MeanFrequency,excelsheet_filename,'Sheet',3,'Range','B2:GG62');%edit if you have more steps/reps
writematrix(ISImin,excelsheet_filename,'Sheet',4,'Range','B2:GG62');%edit if you have more steps/reps
writematrix(MaxIF,excelsheet_filename,'Sheet',5,'Range','B2:GG62');%edit if you have more steps/reps
writematrix(vec_rheobase,excelsheet_filename,'Sheet',6,'Range','B2:GG2');%edit if you have more steps/reps
writematrix(vec_rheobase_frequency,excelsheet_filename,'Sheet',6,'Range','B3:GG3');%edit if you have more steps/reps
writematrix(vec_fmax,excelsheet_filename,'Sheet',7,'Range','B2:GG2');%edit if you have more steps/reps
writematrix(vec_current_fmax,excelsheet_filename,'Sheet',7,'Range','B3:GG3');%edit if you have more steps/reps


