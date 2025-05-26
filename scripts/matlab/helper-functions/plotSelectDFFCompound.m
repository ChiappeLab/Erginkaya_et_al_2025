% This function takes delta F over F (DFF) traces from calcium imaging data
% aligned in CaImagingPopulationAnalysis210718.m (or a more recent .m file)
% and plots specified neurons and optic flow stimuli. Returns figure handle
% as output.

%% Log
% 
% 240130
% Added support to handle latest Trans+Rot stimuli and its controls
% Divided the plots into two, one for controls and one for the compound OF
% 
% 220128
% 
% Added plotType that specifies whether to plot individual traces ('Line')
% or plot the population std as a shaded area ('Shaded')

function f=plotSelectDFFCompound(data,neuronIndex,stimIndex,...
        minFrameRate,colors,prepIndices,prepIndicesFly,gcIndices,...
        gcIndicesFly,sideIndices,sideIndicesFly,Stimuli,regions,limsY,...
        plotType)

cRect = [0 0.5 0 0.1]; % color of the rectangle depicting stimulus motion
lw = 3; % line width of the mean trace
% create a figure with tight_subplot, to minimize empty space in between
f=figure('WindowState','maximized');

switch length(stimIndex)
    % manual reordering based on translational speed
    case 2
        stimOrder = [1,2];
        [ha,~]=tight_subplot(length(stimIndex)/2,length(neuronIndex)*2,...
    0.004,0.025,0.025);
    case 3
        stimOrder = [1,2,3];
        [ha,~]=tight_subplot(3,1,0.004,0.025,0.025);
    case 6
        stimOrder = 1:6;
        [ha,~]=tight_subplot(3,2,0.004,0.025,0.025);
    case 9
        stimOrder = [1,3,2,4,6,5,7,9,8];
        [ha,~]=tight_subplot(3,3,0.004,0.025,0.025);
    case 10
        stimOrder = [5,6,3,4,1,2,7,8,9,10];
        [ha,~]=tight_subplot(length(stimIndex)/2,length(neuronIndex)*2,...
    0.004,0.025,0.025);
    case 12
        stimOrder = [7,8,3,4,1,2,5,6,9,10,11,12];
        [ha,~]=tight_subplot(length(stimIndex)/2,length(neuronIndex)*2,...
    0.004,0.025,0.025);
    case 20
        stimOrder = [1,2,3,4,9,10,13,14,17,18,7,8,5,6,11,12,15,16,19,20];
        [ha,~]=tight_subplot(length(stimIndex)/2,length(neuronIndex)*2,...
    0.004,0.025,0.025);
    case 22
        stimOrder = [1,2,3,4,9,10,15,16,19,20,5,6,11,12,17,18,21,22];
        [ha,~]=tight_subplot(length(stimIndex)/2,length(neuronIndex)*2,...
    0.004,0.025,0.025);
    otherwise
        error('unrecognized number of stimuli, check stimIndex variable')
end

n = length(neuronIndex);
for ii = 1:length(neuronIndex) % for each unique neuron-region
    i = neuronIndex(ii);
    frameNum = size(data{i,1},1);
    TimeAxisPD = linspace(0,frameNum/(minFrameRate/1000),frameNum);
    limsX = [min(TimeAxisPD)+1 ceil(max(TimeAxisPD))-0.5];
    for jj = 1:length(stimIndex) % for each OF stimulus
        j = stimIndex(stimOrder(jj));
        switch length(stimIndex)
            case {3,6,9}
                subIndex = jj;
            otherwise
                subIndex = jj+floor((jj+1)/2)*(2*(n-1))-(2*(n-1))+(2*(ii-1));
        end
        
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
                    nFlies=sum(~isnan(nanmean(data{i,j}(:,iPlot),1)));
                    h=shadedErrorBar(TimeAxisPD,...
                        nanmean(data{i,j}(:,iPlot),2),...
                        nanstd(data{i,j}(:,iPlot),0,2)/...
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
        dashes = strfind(Stimuli{j},'_');
        if length(dashes)==2
            title(Stimuli{j}(1:dashes(2)-1),'Interpreter','none')
        else
            title(Stimuli{j}([1:dashes(1)-1 dashes(2)+1:dashes(3)-1 ...
                dashes(3)+1:end]),'Interpreter','none')
        end
        % add axis labels to the first column, hide axes on other subplots
        if j==stimIndex(stimOrder(1))
            ylabel('Response (Z-score)')
            xlabel('Time (s)')
            if i < length(regions)
               set(gca,'XTick',[])
               ytl=round(linspace(limsY(1),limsY(2),8));
               set(gca,'YTick',ytl)
               set(gca,'YTickLabel',num2cell(ytl,1))
            end
        else
            axis off
        end
        % add the neuron, number of flies (N) and neurites (n) per row
        % if j == stimIndex(stimOrder(2))
        %     ltext = [regions{i} ...
        %         ' n=' num2str(sum(~isnan(nanmean(data{i,j}(:,iPlot),1))))];
        %     th = text(0.1,0.9*limsY(2),ltext);
        %     set(th,'FontSize',16)
        % end
    end
end
