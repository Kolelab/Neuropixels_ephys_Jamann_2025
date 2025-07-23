function tbl = nj_collect_pair_values( db, conditions,pair_fields,fields,crit)
%nj_collect_values. Collects values about pairs from nj_db and returns table
%
%  tbl = nj_collect_pair_values( db, conditions, period_type, fields )
%
%     tbl: Table with subject, condition, val as columns
%
% 2025, Alexander Heimel

if nargin<5 || isempty(crit)
    crit = '';
end
if nargin<4 || isempty(fields)
    fields = {};
end
if nargin<3 || isempty(pair_fields)
    pair_fields = {};
end

if ~iscell(pair_fields) % i.e. single field
    pair_fields = {pair_fields};
end
n_pair_fields = length(pair_fields);

if ~iscell(fields) % i.e. single field
    fields = {fields};
end
n_fields = length(fields);


subject = {};
condition = {};
pair_vals = cell(1,n_pair_fields);
vals = cell(1,n_fields);

count = 1;
for c = 1:2
    ind = find(arrayfun(@(x) strcmp(x.condition,conditions{c}),db));
    for i = 1:length(ind)
        measures = db(ind(i)).measures;
        measures = measures(find_record(measures,crit));
        if  isempty(measures)
            continue
        end
        for m = 1:length(measures)
            pairs = measures(m).pairs;
            if isempty(pairs)
                continue
            end
            n = length(pairs);
            for f = 1:n_pair_fields
                if ischar(pairs(1).(pair_fields{f}))
                    v = {pairs(:).(pair_fields{f})};
                else

                    v = [pairs(:).(pair_fields{f})];  % causes problems if field is empty
                end
                if length(v) ~= n
                    logmsg('Missing values. Solve!')
                end
                pair_vals{f}(count:count+n-1,1) = v(:) ;
            end
            subject(count:count+n-1,1) = {db(ind(i)).subject};
            condition(count:count+n-1,1) = conditions(c);
            for f = 1:n_fields
                vals{f}(count:count+n-1,1) = {measures(m).(fields{f})};
            end
            count = count + n;

            % v = arrayfun(@(p) p.(field),pairs);
            % n = length(v);
            % subject(count:count+n-1,1) = {db(ind(i)).subject};
            % condition(count:count+n-1,1) = conditions(c);
            % val(count:count+n-1,1) = v(:) ;
            % count = count + n;
        end % m
    end % i
end % condition c
subject = categorical(subject);
condition = categorical(condition);
% tbl = table(subject,condition,val);
tbl = table(subject,condition,vals{:},pair_vals{:},'VariableNames',{'subject','condition',fields{:},pair_fields{:}});


