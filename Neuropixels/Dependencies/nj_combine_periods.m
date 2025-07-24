function record = nj_combine_periods(record)
%nj_combine_periods. Combines pairs from different periods
%
% 2025, Alexander Heimel

measures = record.measures;

% remove best period_type
mask = arrayfun(@(m) strcmp(m.period_type,'best') , measures);
measures(mask) = [];
record.measures = measures;


mask = arrayfun(@(m) ~isempty(m.pairs) , measures);


measures = measures(mask);

clusters = unique([measures.cluster]);
n_clusters = length(clusters);

best_measures = measures([]);
count = 1;
for i = 1:n_clusters


    % take opto_test if available
    ind_opto_test = find([measures.cluster]==clusters(i) & contains({measures.period_type},'opto_test'));
    if ~isempty(ind_opto_test)
        best_measures(count) = measures(ind_opto_test);
        best_measures(count).period_type = 'best';
        count = count + 1;
        continue
    end

    % find all measures for this cluster
    ind_measures = find([measures.cluster]==clusters(i));

    ind_measures = find([measures.cluster]==clusters(i) & ...
        (contains({measures.period_type},'whisker_test') | ...
        contains({measures.period_type},'opto_plus_whisker') | ...
        contains({measures.period_type},'all')) )  ;
    if isempty(ind_measures)
        continue
    end

    pairs = [];
    % join all pairs
    for j = ind_measures
        pairs = [pairs measures(j).pairs];
    end
    % get only unique pairs with peak_snr
    new_pairs = [];
    post_clusters = unique([pairs.post_cluster]);
    for j = 1:length(post_clusters)
        ind = find([pairs.post_cluster]==post_clusters(j));

        [~,ind_ind] = max([pairs(ind).peak_snr]);
        new_pairs = [new_pairs pairs(ind(ind_ind))];
    end


    best_measures(count) = measures(ind_measures(1));
    best_measures(count).period_type = 'best';
    best_measures(count).pairs = new_pairs;
    count = count + 1;
end % i

record.measures = [record.measures best_measures];