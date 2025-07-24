function [spont_start,spont_end,spont_duration] = nj_get_spont_interval(record,sAP)
% nj_get_spont_interval. Returns the start and end of the spontaneous activity recordings
%
%   [spont_start,spont_end,spont_duration] = nj_get_spont_interval(record,[sAP])
%
%   spont_start: start of spontaneous period
%   spont_end: end of spontaneous period
%   spont_duration: duration of spontaneous period
%
% Based on getspontrec by Nora Jamann
% 2025, Alexander Heimel

if nargin<2 || isempty(sAP)
    sAP = nj_load_data(record);
end

spont_test = record.spont_test;

if spont_test == -1  % no stimuli were recorded --> entire file
    spont_start = 0;
    spont_end = str2double(sAP.sSources.sMetaNI.fileTimeSecs);
    spont_duration = spont_end - spont_start;
elseif spont_test == 0  %in case of spontaneous recording before first stimulus
    spont_start = min(arrayfun(@(s) s.SpikeTimes(1),sAP.sCluster)); % first spike
    if ~isempty(sAP.cellBlock)
        spont_end = sAP.cellBlock{1, 1}.vecStimOnTime(1);
    else
        spont_end =  max(arrayfun(@(s) s.SpikeTimes(end),sAP.sCluster)); % last spike
    end
    spont_duration = spont_end - spont_start;
elseif mod(spont_test, 1) ~= 0  % not integer, recordings were between two stim recordings
    start_idx = spont_test-0.5;
    end_idx= start_idx + 1;
    spont_start = sAP.cellBlock{1, start_idx}.vecStimOffTime(end)+26;
    spont_end = sAP.cellBlock{1, end_idx}.vecStimOnTime(1);
    spont_duration = spont_end - spont_start;
else % non-zero integer
    spont_start = sAP.cellBlock{1, spont_test}.vecStimOnTime(1) ...
        - sAP.cellBlock{1, spont_test}.dblPrePostWait ...
        - sAP.cellBlock{1, spont_test}.dblPulseWait;
    spont_end = sAP.cellBlock{1, spont_test}.vecStimOffTime(end);
    spont_duration = spont_end - spont_start;
end

