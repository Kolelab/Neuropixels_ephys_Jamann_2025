function [Zeta,ZetaP,Zetapeak] = zeta_baseline(dblBaselineDuration,dblUseMaxDur,vec_StimTimes,vec_SpikeTimes,intPlot)
%Run zeta but with a baseline included

%intPlot = 4;
	%	- intPlot: integer, plotting switch (0=none, 1=inst. rate only, 2=traces only, 3=raster plot as well, 4=adds latencies in raster plot) (default: 0)
intResampNum = 250;
vec_StimTimes_with_baseline = vec_StimTimes - dblBaselineDuration;


[dblZetaP,sZETA,sRate,sLatencies] = zetatest(vec_SpikeTimes,vec_StimTimes_with_baseline,dblUseMaxDur,intResampNum,intPlot);


dblBaselineDurationMs = dblBaselineDuration*1000;
drawnow;hFig = gcf;
for intPlotNr=1:numel(hFig.Children)
	%adjust x-ticks
	if contains(hFig.Children(intPlotNr).XLabel.String,'Time ')
		set(hFig.Children(intPlotNr),'xticklabel',cellfun(@(x) num2str(str2double(x)-dblBaselineDuration),get(hFig.Children(intPlotNr),'xticklabel'),'UniformOutput',false));
	end
	%adjust timings in title
	strTitle = hFig.Children(intPlotNr).Title.String;
	[vecStart,vecStop]=regexp(strTitle,'[=].*?[m][s]');
	for intEntry=1:numel(vecStart)
		strOldNumber=hFig.Children(intPlotNr).Title.String((vecStart(intEntry)+1):(vecStop(intEntry)-2));
		strTitle = strrep(strTitle,strcat('=',strOldNumber,'ms'),strcat('=',num2str(str2double(strOldNumber)-dblBaselineDurationMs),'ms'));
	end
	hFig.Children(intPlotNr).Title.String = strTitle;
end
drawnow;


% sLatencies.Onset = sLatencies.Onset - dblBaselineDuration;
% sLatencies.Peak = sLatencies.Peak - dblBaselineDuration;
% sLatencies.ZETA = sLatencies.ZETA - dblBaselineDuration;
% sLatencies.ZETA_InvSign = sLatencies.ZETA_InvSign - dblBaselineDuration;
% sZETA_pb.vecSpikeT = sZETA_pb.vecSpikeT - dblBaselineDuration;
% sZETA_pb.vecLatencies = sZETA_pb.vecLatencies - dblBaselineDuration;
% sZETA_pb.dblZetaT = sZETA_pb.dblZetaT - dblBaselineDuration;
% sZETA_pb.dblZetaT_InvSign = sZETA_pb.dblZetaT_InvSign - dblBaselineDuration;
% sRate_pb.vecT = sRate_pb.vecT - dblBaselineDuration;
% sRate_pb.dblPeakTime = sRate_pb.dblPeakTime - dblBaselineDuration;
% sRate_pb.dblOnset = sRate_pb.dblOnset - dblBaselineDuration;

Zeta = sZETA.dblZETA;
ZetaP= dblZetaP;
Zetapeak= sZETA.dblPeakT - dblBaselineDuration; 


end