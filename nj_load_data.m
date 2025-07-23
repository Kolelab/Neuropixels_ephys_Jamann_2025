function sAP = nj_load_data(record)
%nj_load_data. Loads Jamann Acquipix sAP data
%
%  sAP = nj_load_data(record)
%
% 2025, Alexander Heimel

params = nj_default_parameters();

filename = fullfile(params.projectfolder,'Data_collection',...
    'NPXdataAlexander',record.reclength,record.condition,[record.sessionid '.mat']);

logmsg(['Loading ' filename]);
load(filename,'sAP');

