function tbl = nj_collect_values( db, conditions, fields, crit )
%nj_collect_values. Collects values from nj_db and returns table
%
%  tbl = nj_collect_values( db, conditions, period_type, fields )
%
%     tbl: Table with subject, condition, val as columns
%
% 2025, Alexander Heimel

if nargin<4 || isempty(crit)
    crit = '';
end

if ~iscell(fields) % i.e. single field
    fields = {fields};
end
n_fields = length(fields);


subject = {};
condition = {};
vals = cell(1,n_fields);

count = 1;
for c = 1:2
    ind = find(arrayfun(@(x) strcmp(x.condition,conditions{c}),db));
    for i = 1:length(ind)
        measures = db(ind(i)).measures;
        measures = measures(find_record(measures,crit));
%        measures = measures(contains({measures.area},area));
        if  isempty(measures)
            continue
        end
        n = length(measures);
        for f = 1:n_fields
            %v = arrayfun(@(m) m.(fields{f}),measures);
            v = [measures(:).(fields{f})];
            if n ~= length(v)
                logmsg('Missing values. Solve');
            end
            vals{f}(count:count+n-1,1) = v(:) ;
        end
        subject(count:count+n-1,1) = {db(ind(i)).subject};
        condition(count:count+n-1,1) = conditions(c);
        count = count + n;
    end % i
end % condition c
subject = categorical(subject);
condition = categorical(condition);
tbl = table(subject,condition,vals{:},'VariableNames',{'subject','condition',fields{:}});

