% This function takes delta F over F (DFF) traces from calcium imaging data
% aligned in CaImagingPopulationAnalysis210718.m (or a more recent .m file)
% and plots specified neurons and optic flow stimuli. Returns figure handle
% as output.

%% Log
% 
% 220128
% 
% Added plotType that specifies whether to plot individual traces ('Line')
% or plot the population std as a shaded area ('Shaded')

function f=plotSelectDS(data1st,data2nd,neuronIndex,stimIndex,...
        minFrameRate,colors,prepIndices,prepIndicesFly,gcIndices,...
        gcIndicesFly,sideIndices,sideIndicesFly,Stimuli,regions,limsY,...
        plotType,normType)

cRect = [0 0.5 0 0.1]; % color of the rectangle depicting stimulus motion
lw = 3; % line width of the mean trace
% create a figure with tight_subplot, to minimize empty space in between
f=figure('WindowState','maximized');
subIndex = 0;
[ha,~]=tight_subplot(length(neuronIndex),length(stimIndex),...
    0.004,0.025,0.025);

for i = neuronIndex(1:length(neuronIndex)) % for each unique neuron-region
    frameNum = size(data1st{i,1},1);
    TimeAxisPD = linspace(0,frameNum/(minFrameRate/1000),frameNum);
%     limsX = [min(TimeAxisPD)+0.5 ceil(max(TimeAxisPD))-0.5];
    limsX = [min(TimeAxisPD)-0.5 ceil(max(TimeAxisPD))-0.5];
    iPlot = prepIndices{i} & gcIndices{i} & sideIndices{i};
    iFly = prepIndicesFly{i} & gcIndicesFly{i} & sideIndicesFly{i};
    for j = stimIndex(1:length(stimIndex)) % for each OF stimulus
        subIndex = subIndex + 1;
        if ~isempty(data1st{i,j}(:,iPlot)) % if there is imaging data
            axes(ha(subIndex)); %#ok<*LAXES>
            % plot individual neurites, and mean across neurites
            if contains(Stimuli{j},{'Lift','YawRightFlipped','Roll','Sideslip'})
                dataToPlot = -data1st{i,j}(:,iPlot)+data2nd{i,j}(:,iPlot);
            else
                dataToPlot = data1st{i,j}(:,iPlot)-data2nd{i,j}(:,iPlot);
            end
            switch plotType
                case 'Line'
                    plot(TimeAxisPD,dataToPlot); hold on;
                    plot(TimeAxisPD,nanmean(dataToPlot,2),...
                        'color',colors(i,:),'LineWidth',lw);
                case 'Shaded'
                    nFlies=sum(~isnan(nanmean(dataToPlot)));
                    h=shadedErrorBar(TimeAxisPD,...
                        nanmean(dataToPlot,2),...
                        nanstd(dataToPlot,0,2)/...
                        sqrt(nFlies));hold on
                    set(h.mainLine,'Color',colors(i,:),...
                        'LineWidth',lw)
                    set(h.patch,'FaceColor',colors(i,:),...
                        'FaceAlpha',0.2)
                    set(h.edge(1),'Color',colors(i,:))
                    set(h.edge(2),'Color',colors(i,:))
                    clearvars h
            end
            rectangle('Position',[2 limsY(1) 2 limsY(2)-limsY(1)],...
                'FaceColor', cRect);
            line(limsX,[0 0],'Color','k','LineWidth',0.5,'LineStyle','--')
            set(gca,'YLim',limsY,'XLim',limsX)
        end
        if i==neuronIndex(1) % add stimulus info to the top of each column
            dashes = strfind(Stimuli{j},'_');
            if length(dashes)==2
                title(Stimuli{j}(4:dashes(1)-1))
            else
                title(Stimuli{j}([4:dashes(1)-1 dashes(2)+1:dashes(3)-1]))
            end
        end
        % add axis labels to the first column, hide axes on other subplots
        if j==stimIndex(1)
            ylabel(['Response (', normType, ')'])
            xlabel('Time (s)')
            if i > neuronIndex(1)
                set(gca,'XTick',[])
            end
            yticklabels('auto')
        else
            axis off
        end
        % add the neuron, number of neurites (n) per row
        if j == stimIndex(ceil(length(stimIndex)/2))
            ltext = [regions{i} ...
                ' n=' num2str(sum(~isnan(nanmean(data1st{i,j}(:,iPlot),1))))];
            th = text(0.1,0.9*limsY(2),ltext);
            set(th,'FontSize',16)
        end
    end
end
