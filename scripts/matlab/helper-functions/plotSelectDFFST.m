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

function f=plotSelectDFFST(data,neuronIndex,stimIndex,...
        minFrameRate,colors,prepIndices,prepIndicesFly,gcIndices,...
        gcIndicesFly,sideIndices,sideIndicesFly,Stimuli,regions,limsY,...
        plotType)

cRect = [0 0.5 0 0.1]; % color of the rectangle depicting stimulus motion
lw = 3; % line width of the mean trace
% create a figure with tight_subplot, to minimize empty space in between
f=figure('WindowState','maximized');
subIndex = 0;
[ha,~]=tight_subplot(length(neuronIndex),length(stimIndex),...
    0.004,0.025,0.025);

for ii = 1:length(neuronIndex) % for each unique neuron-region
    i = neuronIndex(ii);
    frameNum = size(data{i,1},1);
    TimeAxisPD = linspace(0,frameNum/(minFrameRate/1000),frameNum);
    limsX = [min(TimeAxisPD)+1 ceil(max(TimeAxisPD))];
    for jj = 1:length(stimIndex) % for each OF stimulus
        switch regions{i}
            case 'H2axon'
                j = stimIndex(jj)+8;
                ji = stimIndex(ceil(length(stimIndex)/2))+8;
            otherwise
                j = stimIndex(jj);
                ji = stimIndex(ceil(length(stimIndex)/2));
        end
        
        subIndex = subIndex + 1;
        iPlot = prepIndices{i} & gcIndices{i} & sideIndices{i};
        iPlot = true(size(data{i,j},2),1);
        iFly = prepIndicesFly{i} & gcIndicesFly{i} & sideIndicesFly{i};
        iFly=iPlot;
        if ~isempty(data{i,j}(:,iPlot)) % if there is imaging data
            axes(ha(subIndex)); %#ok<*LAXES>
            % plot individual neurites, and mean across neurites, 1st epoch
            switch plotType
                case 'Line'
                    plot(TimeAxisPD,data{i,j}(:,iPlot)); hold on;
                    plot(TimeAxisPD,nanmean(data{i,j}(:,iPlot),2),...
                        'color',colors(i,:),'LineWidth',lw);
                case 'Shaded'
                    h=shadedErrorBar(TimeAxisPD,...
                        nanmean(data{i,j}(:,iPlot),2),...
                        nanstd(data{i,j}(:,iPlot),0,2)/...
                        sqrt(sum(iPlot)));hold on
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
                title(Stimuli{j}([4:dashes(1)-1 dashes(2)+1:dashes(3)-1 ...
                    dashes(3)+1:end]))
            end
        end
        % add axis labels to the first column, hide axes on other subplots
        if j==stimIndex(1)
            ylabel('Response')
            xlabel('Time (s)')
            if i < length(regions)
                set(gca,'XTick',[])
            end
        else
            axis off
        end
        % add the neuron, number of flies (N) and neurites (n) per row
        
        if j == ji
            ltext = [regions{i} ...
                ' n=' num2str(sum(~isnan(nanmean(data{i,j}(:,iPlot),1))))];
            th = text(0.1,0.9*limsY(2),ltext);
            set(th,'FontSize',16)
        end
    end
end
