% This function takes traces from calcium imaging data, calculates and
% plots the peak response for specified neurons and optic flow stimuli.
% Returns figure handle as output.

function f=plotSelectPeak(data1st,data2nd,neuronIndex,stimIndex,...
        colors,prepIndices,prepIndicesFly,gcIndices,gcIndicesFly,...
        sideIndices,sideIndicesFly,Stimuli,regions,plotType,...
        groupType,iPeak)      

lw = 3; % line width of the mean trace
sp = 0.5; % separation between each box plot
switch groupType
    case 'neuron'
        h = cell(length(neuronIndex),length(stimIndex));
        % initialize x-axis matrix
        x = zeros(length(neuronIndex),length(stimIndex)*2);
        % initialize maximum y-axis value
        my = zeros(1,length(stimIndex));
        % initialize legend handle and text
        lh = gobjects(1,length(stimIndex));
        ltext = cell(1,length(stimIndex));
    case 'stimulus'
        h = cell(length(neuronIndex),length(stimIndex));
        % initialize x-axis matrix
        x = zeros(length(neuronIndex),length(stimIndex)*2);
        % initialize maximum y-axis value
        my = zeros(1,length(neuronIndex));
        % initialize legend handle and text
        lh = gobjects(1,length(neuronIndex));
        ltext = cell(1,length(neuronIndex));
    otherwise
        error(['Grouping not specified correctly. ',...
            'Use groupType = [''neuron'' or ''stimulus'']'])
end

% generate the figure
f=figure('WindowState','maximized');
for i = 1:length(neuronIndex) % for each unique neuron-region
    ni = neuronIndex(i); % nauron index in the original data matrix
    iPlot = prepIndices{ni} & gcIndices{ni} & sideIndices{ni};
    iFly = prepIndicesFly{ni} & gcIndicesFly{ni} & sideIndicesFly{ni};
    for j = 1:length(stimIndex) % for each OF stimulus
        si = stimIndex(j); % stimulus index in the original data matrix
        switch groupType
            case 'neuron'
                x(i,j) = (i-1)*length(stimIndex) + 1 + (j-1)*sp;
                x(i,j+length(stimIndex)) = length(stimIndex)*...
                    i + 1 + (j-1)*sp;
            case 'stimulus'
                x(i,j) = (j-1)*length(neuronIndex) + 1 + (i-1)*sp;
                x(i,j+length(stimIndex)) = (j+length(stimIndex)-1)*...
                    length(neuronIndex) + 1 + (i-1)*sp;
        end
        if ~isempty(data1st{ni,si}(:,iPlot)) % if there is data
            dataToPlot1 = mean(data1st{ni,si}(iPeak,iPlot));
            switch plotType
                % Plot each fly (dots), mean (line) and std (box)
                case 'Box'
                    h{i,j} = notBoxPlot(dataToPlot1,x(i,j));
                    set(h{i,j}.data,'MarkerFaceColor',colors(ni,:))
            end
            hold on
        end
        if ~isempty(data2nd{ni,si}(:,iPlot)) % if there is data
            dataToPlot2 = mean(data2nd{ni,si}(iPeak,iPlot));
            switch plotType
                % Plot each fly (dots), mean (line) and std (box)
                case 'Box'
                    h{i,j} = notBoxPlot(dataToPlot2,...
                        x(i,j+length(stimIndex)));
                    set(h{i,j}.data,'MarkerFaceColor',colors(ni,:))
            end
            hold on
        end
        % save the maximum y values for significance plot alignment
        maxy = max(nanmax(dataToPlot1),nanmax(dataToPlot2));
        if isempty(maxy)
            maxy = 0;
        end
        my(i) = max(my(i),maxy);
        % add the neuron, number of flies (i) and neurites (i) for legend
        lh(1,i) = plot(nan,nan,'-','Color',colors(ni,:),'LineWidth',lw);
        ltext{1,i} = [regions{ni} ' N=' num2str(sum(iFly)) ...
            ' n=' num2str(sum(iPlot))];
    end
end


limsX = xlim;
line(limsX,[0 0],'Color','k','LineWidth',0.5,'LineStyle','--')
ylabel('response (peak)')

switch groupType
    case 'neuron'
        % locations for x-tick marks
        iXtick = sort(reshape(x,1,numel(x)));
        set(gca,'XTick',iXtick,'XTickLabel',regions(neuronIndex),...
            'FontSize',12,'FontName','Arial',...
            'XTickLabelRotation',45,'TickLabelInterpreter','none')
    case 'stimulus'
        % locations for x-tick marks
        iXtick = x(round(size(x,1)/2),:);
        set(gca,'XTick',iXtick,'XTickLabel',Stimuli(stimIndex),...
            'FontSize',12,'FontName','Arial',...
            'XTickLabelRotation',45,'TickLabelInterpreter','none')
end

for j = 1:length(stimIndex) % for each OF stimulus
    si = stimIndex(j); % stimulus index in the original data matrix
    % Perform Mann-Whitney-Wilcoxon non parametric test and plot
    % significance indicators if p<0.05
    nck = nchoosek(1:length(neuronIndex),2);
    nph = size(nck,1); % number of comparisons made (used for adjusting the
    % minimum P value required for significance - a.k.a. Bonferroni correction)
    %This is all done just to help with the plots
    nck = [nck; nck(end,:)];
    sameComparison = 1;
    for ii = 1 : size(nck,1) - 1
        if nck(ii+1) == nck(ii)
            sameComparison = sameComparison + 1;
        else
            sameComparison = 1;
        end
    end
    for ii = 1 : size(nck,1) - 1
        ni1 = neuronIndex(nck(ii,1));
        ni2 = neuronIndex(nck(ii,2));
        iPlot1 = prepIndices{ni1} & gcIndices{ni1} & sideIndices{ni1};
        iPlot2 = prepIndices{ni2} & gcIndices{ni2} & sideIndices{ni2};
        d1 = mean(data1st{ni1,si}(iPeak,iPlot1)); d1(isnan(d1)) = [];
        d2 = mean(data1st{ni2,si}(iPeak,iPlot2)); d2(isnan(d2)) = [];
        mym = max(my(nck(ii,1)),my(nck(ii,2)));
        stats{ii} = mwwtest(d1,d2,0);
        if stats{ii}.p < 0.05/nph
            % uncomment below if you want to display the p-values per pair
            disp([Stimuli{si},'_1-',regions{ni1},'-',regions{ni2},'-',...
                num2str(stats{ii}.p(2))])
            sigline([x(nck(ii,1),j),x(nck(ii,2),j)],[],mym+mym*0.01*ii);
        end
        d1 = mean(data2nd{ni1,si}(iPeak,iPlot1)); d1(isnan(d1)) = [];
        d2 = mean(data2nd{ni2,si}(iPeak,iPlot2)); d2(isnan(d2)) = [];
        mym = max(my(nck(ii,1)),my(nck(ii,2)));
        stats{ii} = mwwtest(d1,d2,0);
        if stats{ii}.p < 0.05/nph
            % uncomment below if you want to display the p-values per pair
            disp([Stimuli{si},'_2-',regions{ni1},'-',regions{ni2},'-',...
                num2str(stats{ii}.p(2))])
            sigline([x(nck(ii,1),j+length(stimIndex)),x(nck(ii,2),...
                j+length(stimIndex))],[],mym+mym*0.01*ii);
        end
    end
end
legend(lh,ltext,'Location','Best')
