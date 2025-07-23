function tbl = nj_all_data_comparison(tbl,field,label,scale,clp,yl,show_points,handle)
%nj_all_data_comparison. Plots comparison of Ctrl and Cpz for all data
%
% 2025, Alexander Heimel

if nargin<8 || isempty(handle)
    handle = [];
end

if nargin<4 || isempty(scale)
    scale = 'linear';
end
if nargin<5 || isempty(clp)
    clp = [];
end
if nargin<6 || isempty(yl)
    yl = [];
end
if nargin<7 || isempty(show_points)
    show_points = 1;
end


tbl(isnan(tbl.(field)),:) = [];

if ~isempty(clp)
    tbl.(field)(tbl.(field)<clp(1)) = clp(1);
    tbl.(field)(tbl.(field)>clp(2)) = clp(2);
end

mean_field = ['mean_' field];
mean_mean_field = ['mean_mean_' field];

tbl_per_subject = groupsummary(tbl,{'subject', 'condition'},"mean",field);
tbl_per_condition = groupsummary(tbl_per_subject,'condition',"mean",mean_field);

switch scale
    case 'log'
        tbl.(field) = log10(tbl.(field));
        tbl_per_subject.(mean_field) = log10(tbl_per_subject.(mean_field));
        tbl_per_condition.(mean_mean_field) = log10(tbl_per_condition.(mean_mean_field));
end

if isnumeric(tbl.(field))
    lme = fitlme(tbl,[field ' ~ condition + (1|subject)']);
    p_value = coefTest(lme);
    test = 'Nested t-test';
else
    logmsg([field ' is not numeric. Skipping LME and doing two-sample test on means per subject instead.'] );
    [~,p_value] = ttest2(tbl_per_subject.(mean_field)(tbl_per_subject.condition=="Ctrl"),...
        tbl_per_subject.(mean_field)(tbl_per_subject.condition=="Cpz") );
    test = 'Unpaired t-test';
end

% all units
if show_points
    
    ivt_graph({tbl.(field)(tbl.condition=="Ctrl"),tbl.(field)(tbl.condition=="Cpz")},[],...
        'style','level','spaced',3,'test','none','markersize',3,'barwidth',0,...
        'errorbars','none','markers','open_circle','color',[1 0.7 0.7],'axishandle',handle);
    c = get(gca,'Children');
    for i=1:length(c)
        set(c(i),'Color',[0.7 0.7 0.7])
    end
else
    figure
end


% all mice
ivt_graph({...
    tbl_per_subject.(mean_field)(tbl_per_subject.condition=="Ctrl"),...
    tbl_per_subject.(mean_field)(tbl_per_subject.condition=="Cpz")},[],...
    'color',{[1 1 1],[0 0 0]},'test','none','errorbars','none',...
    'style','level','spaced',3,'axishandle',gca,'markersize',6,'barwidth',0);

% all conditions
ivt_graph({...
    tbl_per_condition.(mean_mean_field)(tbl_per_condition.condition=="Ctrl"),...
    tbl_per_condition.(mean_mean_field)(tbl_per_condition.condition=="Cpz")},[],...
    'color',{[1 1 1],[0 0 0]},'test','none',...
    'style','level','spaced',3,'axishandle',gca,'markersize',6,...
    'barwidth',0.5,'showpoints',0);

fontsize(scale=0.8);
ylabel(label)

switch scale
    case 'log'
        set(gca,'ytick',[-2 -1 0 1 2],'yticklabel',[0.01 0.1 1 10 100]);
        ylim([-4 2])
end

set(gca,'xtick',[1 2],'xticklabel',{'Ctrl','Cpz'});

set(gca,'units','centimeters','position',[3 2 1.5 3.5]);

if ~isempty(yl)
    ylim(yl);
end

yl = ylim();
y = yl(1) + 0.9 *(yl(2)-yl(1));
plot([1 2],y*[1 1],'k-','LineWidth',1)
star = 'ns';
if p_value<0.001
    star = '***';
elseif p_value<0.01
    star = '**';
elseif p_value<0.05
    star = '*';
end
text(1.5,y,star,'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',10)



disp([test ' P = ' num2str(p_value,2) ...
    '. n = ' num2str(sum(tbl.condition=="Ctrl")) ...
    ' neurons/pairs, N = ' num2str(sum(tbl_per_subject.condition=="Ctrl")) ' mice (Ctrl)' ...
    ', n = ' num2str(sum(tbl.condition=="Cpz")) ...
    ' neurons/pairs, N = ' num2str(sum(tbl_per_subject.condition=="Cpz")) ...
    ' mice (Cpz).']);
