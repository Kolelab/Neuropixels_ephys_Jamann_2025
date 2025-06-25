function [matFrequencies_area] = frequenciesperarea(area_names,AreaIdxList,frequencies)

%Inputs
%area_names = vector with names of all areas 
%AreaIdxList = vector with indices of areas that cluster belongs to
%frequencies = mean frequencies of all clusters
%Outputs
%Frequencies array = array of all mean frequencies of good clusters divided
%by area

%%
%Generate a table with all individual frequencies of each
%cluster,divided by area, for all animals
%Find frequencies and waveform parameters of all clusters of a certain area
%generate matrix with all areas

num_areas = numel(area_names);
num_clust_area = zeros(1,num_areas);

% find number of clusters per area to determine the array size
for areaIdx = 1:num_areas
    num_clust_area(areaIdx) = numel(find(AreaIdxList == areaIdx));
end

%preallocate 
max_num = max(num_clust_area);
matFrequencies_area = nan(max_num, num_areas);

%index from full list of clusters, divide by different areas
for area_idx = 1:num_areas    
  Idx_clust_in_area = find(AreaIdxList == area_idx);
  matFrequencies_area(1:num_clust_area(area_idx),area_idx) = frequencies(Idx_clust_in_area)';
end
end

