% This script takes pre-processed behavior data and generates the plots
% shown in the manuscript

clearvars

%% User Input

% Get the full path of this script
scriptFullPath = mfilename('fullpath');
% Directory of the script: repository/scripts/matlab
scriptDir = fileparts(scriptFullPath);

% Go up two levels to reach the repository root:
%   fileparts(scriptDir) gives repository/scripts
%   fileparts(fileparts(scriptDir)) gives repository
repoDir = fileparts(fileparts(scriptDir));

% Build the data directory and results directory paths relative to the repository root
baseDir = fullfile(repoDir, 'data', 'Behavior\');


% Option flag to save figures (1 to save, 0 to just display)
SaveFigures = 0;

% plotting order
gOrder = {'None','Both','Left','Right'};

% plotting colors
cReg = [0.5273 0.5195 0; 0 0.4 0]; % olive (DNp15) ; green (bIPS)

% use bouts that occur only at the center of the arena
params = GetParamsArenaCenter(); % load parameters and thresholds

%% Figure 7b and 7d

% DNp15 hsFlp Kir subfolder
behDir = 'DNp15Kir\';
% experimental condition
cond = 'Light';
genotypes = {'No cells','Single cell Left','Single Cell Right','Double cell'};
gTypes = {'None','Left','Right','Both'}; % short names for plot labels
saveDir = [baseDir,'Analysis\DNp15Kir\'];

[~,gi] = ismember(gOrder,gTypes); % plotting order indices


% load angular deviation data. If AD.mat doesn't exist, calculate angular
% deviations, store in AD.mat and save it on disk
if ~exist([saveDir,'AngDevVars.mat'],'file')
    % initialize the structures. Each structure will have as many columns
    % as there are genotypes. Cell arrays within each coulm will have size
    % of MxN where M is the number of flies and N is number of total
    % conditions (dark + 8 different light conditions as of 220131).
    % Sequence of light conditions will be stored in the seq variable below
    lightCondNum = 2;
    VfwThrS = 14; % mm/s | the maximum forward velocity that will be
    % considered as slow walking
    VfwThrF = 18; % mm/s | the minimum forward velocity that will be
    % considered as fast walking
    AD = struct();
    AD.VfwThr = [VfwThrS,VfwThrF];
    % structure fields for all forward bouts
    AD.Wmu = cell(1,length(genotypes)); % weighted mean
    AD.Wstd = cell(1,length(genotypes)); % std of the weighted mean
    AD.AbsWmu = cell(1,length(genotypes)); % weighted mean of abs angdev
    AD.BWmu = cell(1,length(genotypes)); %weighted mean of binarized angdev
    AD.BlSum = cell(1,length(genotypes)); % bout length summed per fly
    AD.Bl = cell(1,length(genotypes)); % bout length per bout per fly
    AD.Adev = cell(1,length(genotypes)); % bout length per bout per fly
    AD.Blmu = cell(1,length(genotypes)); % mean bout length per fly
    AD.Vfmu = cell(1,length(genotypes)); % mean Vfw per bout per fly
    AD.VFperGen = cell(1,length(genotypes)); % Vfw values per genotype
    % structure fields for "fast" forward bouts
    AD.WmuF = cell(1,length(genotypes)); % weighted mean
    AD.WstdF = cell(1,length(genotypes)); % std of the weighted mean
    AD.AbsWmuF = cell(1,length(genotypes)); % weighted mean of abs angdev
    AD.BlSumF = cell(1,length(genotypes)); % bout length summed per fly
    AD.BlF = cell(1,length(genotypes)); % bout length per bout per fly
    AD.BlmuF = cell(1,length(genotypes)); % mean bout length per fly
    AD.VfmuF = cell(1,length(genotypes)); % mean Vfw per bout per fly
    % structure fields for "slow" forward bouts
    AD.WmuS = cell(1,length(genotypes)); % weighted mean
    AD.WstdS = cell(1,length(genotypes)); % std of the weighted mean
    AD.AbsWmuS = cell(1,length(genotypes)); % weighted mean of abs angdev
    AD.BlSumS = cell(1,length(genotypes)); % bout length summed per fly
    AD.BlS = cell(1,length(genotypes)); % bout length per bout per fly
    AD.BlmuS = cell(1,length(genotypes)); % mean bout length per fly
    AD.VfmuS = cell(1,length(genotypes)); % mean Vfw per bout per fly
    % All angular deviation values concatenated per genotype
    AD.AdevAll = cell(1,length(genotypes));

    for g = 1:length(genotypes) % for each genotype

        % directory where the files are located
        path = [baseDir,behDir,cond,'\',genotypes{g},'\'];
        disp([behDir,cond,'\',genotypes{g},'\'])
        % load the low resolution data relative to the forward segments
        % together with the stimuli sequence
        [VForwSeg,seq] = GetForwSegVData(path, params);
        % NOTE: Fw speed cutoff is stored in params.vft and it can be
        % changed based on slow/fast forward speed analysis preference

        % binning windows for the angular displacements (deg/mm)
        DevCents = [-1000 -3.75 -2.5 -1.25 0 1.25 2.5 3.75 1000];

        % initialize variables
        VR = cell(length(DevCents)-1,1); % rotational velocity
        VF = cell(length(DevCents)-1,1); % forward velocity
        VS = cell(length(DevCents)-1,1); % side velocity
        NB = cell(length(DevCents)-1,1); % number of bouts
        angDevAll = [];
        % Length of each forward walking bout per fly
        Blength = cell(size(VForwSeg,2),size(VForwSeg,1));
        % Angular deviation per bout per fly
        Adev = cell(size(VForwSeg,2),size(VForwSeg,1)); % signed
        AdevAbs = cell(size(VForwSeg,2),size(VForwSeg,1)); % absolute
        AdevB = cell(size(VForwSeg,2),size(VForwSeg,1)); % binarized
        % Mean forward velocity per bout per fly
        Vfw = cell(size(VForwSeg,2),size(VForwSeg,1));

        for jk = 1:length(DevCents)-1 % for each angular deviation bin
            % define boundaries of each bin
            minAngDev = DevCents(jk);
            maxAngDev = DevCents(jk+1);

            for s = 1:size(VForwSeg,1) % for each stimulus
                for n = 1:size(VForwSeg,2) % for each fly
                    for i = 1:length(VForwSeg{s,n}) % for each fw bout
                        if ~isempty(VForwSeg{s,n}{i})% if there is data
                            % low res speed data for this fw segment
                            vr = VForwSeg{s,n}{i}.VrLR; %Vrotational
                            vf = VForwSeg{s,n}{i}.VfLR; %Vforward
                            vs = VForwSeg{s,n}{i}.VsLR; %Vsideways
                            vt = sqrt(vf.*vf + vs.*vs); %Vtranslational
                            angDev = mean(vr./vt); % angular deviation
                            angDevB = angDev/abs(angDev); % binarized

                            %if the angular deviation belong to the bin
                            if angDev > minAngDev && angDev < maxAngDev
                                % store the data appropriately
                                VR{jk} = vertcat(VR{jk}, vr);
                                VF{jk} = vertcat(VF{jk}, vf);
                                VS{jk} = vertcat(VS{jk}, vs);
                                NB{jk} = vertcat(NB{jk}, length(vs));
                                angDevAll = vertcat(angDevAll,angDev);
                                Adev{n,s} = vertcat(Adev{n,s},angDev);
                                AdevB{n,s} = vertcat(AdevB{n,s},angDevB);
                                AdevAbs{n,s} = vertcat(AdevAbs{n,s},...
                                    abs(angDev));
                                Vfw{n,s} = vertcat(Vfw{n,s},mean(vf));
                                Blength{n,s} = vertcat(Blength{n,s},...
                                    length(vs));
                            end
                        end
                    end
                end
            end
        end

        % initialize variables
        % Angular deviation per fly weighted by ALL fw bout length
        AngDevWmu = nan(size(Adev));
        AngDevWvar = nan(size(Adev));
        AngDevAbsWmu = nan(size(Adev));
        AngDevBWmu = nan(size(AdevB));
        BlSum = nan(size(Adev));
        BlW = nan(size(Adev));
        AD.AdevAll{g} = vertcat(AD.AdevAll{g},angDevAll);
        % Angular deviation per fly weighted by FAST fw bout length
        AngDevWmuF = nan(size(Adev));
        AngDevWvarF = nan(size(Adev));
        AngDevAbsWmuF = nan(size(Adev));
        BlSumF = nan(size(Adev));
        AdevF = cell(size(Adev));
        AdevAbsF = cell(size(AdevAbs));
        BlengthF = cell(size(Blength));
        % Angular deviation per fly weighted by SLOW fw bout length
        AngDevWmuS = nan(size(Adev));
        AngDevWvarS = nan(size(Adev));
        AngDevAbsWmuS = nan(size(Adev));
        BlSumS = nan(size(Adev));
        AdevS = cell(size(Adev));
        BlengthS = cell(size(Blength));
        AdevAbsS = cell(size(AdevAbs));

        % logical index arrays for slow and fast bouts
        iFast =cellfun(@(x) x >= VfwThrF, Vfw, 'UniformOutput', false);
        iSlow =cellfun(@(x) x < VfwThrS, Vfw, 'UniformOutput', false);

        for n = 1:size(Adev,1) % for each fly
            for i = 1:size(Adev,2) % for each stimulus

                % Angular deviations of all bouts calculated as a
                % weighted mean based on bout lengths
                AngDevWmu(n,i)=...
                    nansum((Adev{n,i}.*Blength{n,i})/...
                    nansum(Blength{n,i}));
                AngDevAbsWmu(n,i)=...
                    nansum((AdevAbs{n,i}.*Blength{n,i})/...
                    nansum(Blength{n,i}));
                AngDevBWmu(n,i)=...
                    nansum((AdevB{n,i}.*Blength{n,i})/...
                    nansum(Blength{n,i}));
                AngDevWvar(n,i)=...
                    nansum((((Adev{n,i}-AngDevWmu(n,i)).^2)...
                    .*Blength{n,i})/nansum(Blength{n,i}));
                BlSum(n,i) = nansum(Blength{n,i});

                % Angular deviations of fast bouts calculated as a
                % weighted mean based on bout lengths
                AdevF{n,i} = Adev{n,i}(iFast{n,i}); % fast bouts only
                AdevAbsF{n,i} = AdevAbs{n,i}(iFast{n,i});
                BlengthF{n,i} = Blength{n,i}(iFast{n,i});
                AngDevWmuF(n,i)=...
                    nansum((AdevF{n,i}.*BlengthF{n,i})/...
                    nansum(BlengthF{n,i}));
                AngDevAbsWmuF(n,i)=...
                    nansum((AdevAbsF{n,i}.*BlengthF{n,i})/...
                    nansum(BlengthF{n,i}));
                AngDevWvarF(n,i)=...
                    nansum((((AdevF{n,i}-AngDevWmuF(n,i)).^2)...
                    .*BlengthF{n,i})/nansum(BlengthF{n,i}));
                BlSumF(n,i) = nansum(BlengthF{n,i});

                % Angular deviations of slow bouts calculated as a
                % weighted mean based on bout lengths
                AdevS{n,i} = Adev{n,i}(iSlow{n,i}); % slow bouts only
                AdevAbsS{n,i} = AdevAbs{n,i}(iSlow{n,i});
                BlengthS{n,i} = Blength{n,i}(iSlow{n,i});
                AngDevWmuS(n,i)=...
                    nansum((AdevS{n,i}.*BlengthS{n,i})/...
                    nansum(BlengthS{n,i}));
                AngDevAbsWmuS(n,i)=...
                    nansum((AdevAbsS{n,i}.*BlengthS{n,i})/...
                    nansum(BlengthS{n,i}));
                AngDevWvarS(n,i)=...
                    nansum((((AdevS{n,i}-AngDevWmuS(n,i)).^2)...
                    .*BlengthS{n,i})/nansum(BlengthS{n,i}));
                BlSumS(n,i) = nansum(BlengthS{n,i});
            end
        end

        % remove zeros from the data
        AngDevWmu(BlSum==0) = nan;
        AngDevWvar(BlSum==0) = nan;
        AngDevAbsWmu(BlSum==0) = nan;
        AngDevBWmu(BlSum==0) = nan;
        BlSum(BlSum==0) = nan;
        AngDevWmuF(BlSumF==0) = nan;
        AngDevWvarF(BlSumF==0) = nan;
        AngDevAbsWmuF(BlSumF==0) = nan;
        BlSumF(BlSumF==0) = nan;
        AngDevWmuS(BlSumS==0) = nan;
        AngDevWvarS(BlSumS==0) = nan;
        AngDevAbsWmuS(BlSumS==0) = nan;
        BlSumS(BlSumS==0) = nan;
        % standard deviation of the weighted means
        AngDevWstd = sqrt(AngDevWvar);
        AngDevWstdF = sqrt(AngDevWvarF);
        AngDevWstdS = sqrt(AngDevWvarS);

        % horzcat concatenates in the opposite of what we want when the
        % initial array is empty, so we swap the positions of variables
        % during the first concatenation (which happens in 'Dark'
        % condition). If you get errors here, make sure that the first
        % element of conds variable is 'Dark'.
        % Only process the 'Light' condition for this section
        % The switch(cond{c}) is not needed since only 'Light' is used here.
        % The following block is equivalent to the 'Light' case:

        % Commented out the switch/case above for clarity.
        AD.Wmu{g} = horzcat(AD.Wmu{g},AngDevWmu);
        AD.Wstd{g} = horzcat(AD.Wstd{g},AngDevWstd);
        AD.AbsWmu{g} = horzcat(AD.AbsWmu{g},AngDevAbsWmu);
        AD.BWmu{g} = horzcat(AD.BWmu{g},AngDevBWmu);
        AD.BlSum{g} = horzcat(AD.BlSum{g},BlSum);
        AD.Bl{g} = horzcat(AD.Bl{g},Blength);
        AD.Adev{g} = horzcat(AD.Adev{g},Adev);
        AD.WmuF{g} = horzcat(AD.WmuF{g},AngDevWmuF);
        AD.WstdF{g} = horzcat(AD.WstdF{g},AngDevWstdF);
        AD.AbsWmuF{g} = horzcat(AD.AbsWmuF{g},AngDevAbsWmuF);
        AD.BlSumF{g} = horzcat(AD.BlSumF{g},BlSumF);
        AD.BlF{g} = horzcat(AD.BlF{g},BlengthF);
        AD.WmuS{g} = horzcat(AD.WmuS{g},AngDevWmuS);
        AD.WstdS{g} = horzcat(AD.WstdS{g},AngDevWstdS);
        AD.AbsWmuS{g} = horzcat(AD.AbsWmuS{g},AngDevAbsWmuS);
        AD.BlSumS{g} = horzcat(AD.BlSumS{g},BlSumS);
        AD.BlS{g} = horzcat(AD.BlS{g},BlengthS);
        AD.seq = seq'; % add stimulus info (no 'Dark' for Light-only analysis)
        
        %% Figure 7b - Plot example trajectories
        xT = cell(1,1);
        yT = cell(1,1);
        indd = cell(1,1);
        a = 1;
        cmap1 = jet(length(VR));
        xTall = [];
        yTall = [];
        % iterate through the angular deviation bins
        for i = 2 : length(NB)-1
            nbs = NB{i};
            nbScs = [1; cumsum(nbs)];
            % iterate through the forward segments
            for j = 1 : length(nbs)
                vr = VR{i}(nbScs(j):nbScs(j+1))/60;
                vf = VF{i}(nbScs(j):nbScs(j+1))/60;
                vs = VS{i}(nbScs(j):nbScs(j+1))/60;
                x = [];
                y = [];
                th = [];
                if length(vr) < 180
                    % get centered X and Y
                    for ij = 1 : length(vr)
                        if isempty(x)
                            th = vertcat(th,pi*vr(ij)/180);
                            x = vertcat(x,vf(ij)*sin(0)+vs(ij)*cos(0));
                            y = vertcat(y,vf(ij)*cos(0)-vs(ij)*sin(0));
                        else
                            th = vertcat(th,th(end)+pi*vr(ij)/180);
                            x = vertcat(x,x(end)+...
                                vf(ij)*sin(th(ij-1))...
                                +vs(ij)*cos(th(ij-1)));
                            y = vertcat(y,y(end)+...
                                vf(ij)*cos(th(ij-1))...
                                -vs(ij)*sin(th(ij-1)));
                        end
                    end
                end
                xT{a} = x;
                yT{a} = y;
                indd{a} = i;
                if a == 1
                    xTall = x';
                    yTall = y';
                else
                    if length(x) > size(xTall,2)
                        xTall(:,(size(xTall,2)+1):length(x)) = NaN;
                        yTall(:,(size(yTall,2)+1):length(y)) = NaN;
                        xTall(a,:) = x';
                        yTall(a,:) = y';
                    else
                        xTall(a,:) = NaN;
                        yTall(a,:) = NaN;
                        xTall(a,1:length(x)) = x';
                        yTall(a,1:length(y)) = y';
                    end
                end
                a = a + 1;
            end
        end
        yTallBin = 0:round(max(nanmax(yTall)));
        xTallBin = nan(1,length(yTallBin));
        yTallR = round(yTall);
        for i = 1:(round(max(nanmax(yTall)))+1)
            xTallBin(i)=nanmean(xTall(yTallR==i-1));
        end
        % plot a color coded sample of the forward segments
        inds = randperm(length(xT));
        frac = 100; % percentage of data to be plotted
        ythr = 10; % threshold of forward translation (mm in y axis)
        inds = inds(1:round(length(xT)*frac/100));

        fs = 18; % font size
        f=figure('WindowState','maximized');
        hold on
        for i = 1 : length(inds)
            if max(yT{inds(i)})>ythr
                plot(yT{inds(i)}, xT{inds(i)},...
                    'color', cmap1(indd{inds(i)},:),'LineWidth',1.5)
            end
        end
        plot(yTallBin,xTallBin,'color','k','LineWidth',2)

        axis([-2 30 -25 25])
        limsX = xlim;
        line(limsX,[0 0],'Color','k','LineWidth',0.5,'LineStyle','--')
        %             title(['bIPS-',gTypes{g}]) % add title, tick and axis labels
        %             title(['DNp15-',genotypes{g}]) % add title, tick and axis labels
        title(['DNp15-RicinA-',gTypes{g}]) % add title, tick and axis labels
        ylabel('X position (mm)'); xlabel('Y position (mm)')
        set(gca,'FontSize',fs)
        caxis([DevCents(2) DevCents(end-1)])
        colormap(cmap1)
        ch = colorbar('Ticks',DevCents(2:end-1));
        ch.Label.String = 'Ang. Dev. (deg/mm)';
        ch.Label.Rotation = 270;
        ts = timestamp(); % get current date
        if SaveFigures
            set(gcf,'renderer','Painters')
            print(f,strcat(saveDir,ts,'-Paths-',...
                genotypes{g},'-',cond{c},'-',num2str(frac),...
                'perc-yThr',num2str(ythr)),'-dtiff') %#ok<*UNRCH>
            print(f,strcat(saveDir,ts,'-Paths-',...
                genotypes{g},'-',cond{c},'-',num2str(frac),...
                'perc-yThr',num2str(ythr)),'-depsc')
            close(f)
        end
    end

    %     AD.AdevAll = AdevAll;
    for g = 1:length(genotypes) % for each genotype
        nFly = size(AD.Bl{1,g}(:,1),1);
        Blmean = nan(nFly,size(AD.seq,2));
        BlmeanF = nan(nFly,size(AD.seq,2));
        BlmeanS = nan(nFly,size(AD.seq,2));
        for i = 1:size(AD.seq,2) % for each condition
            for n = 1:nFly % for each fly
                Blmean(n,i) = nanmean(AD.Bl{1,g}{n,i});
                BlmeanF(n,i) = nanmean(AD.BlF{1,g}{n,i});
                BlmeanS(n,i) = nanmean(AD.BlS{1,g}{n,i});
            end
        end
        AD.Blmu{1,g} = Blmean;
        AD.BlmuF{1,g} = BlmeanF;
        AD.BlmuS{1,g} = BlmeanS;
    end
    if ~exist(saveDir, 'dir')
        mkdir(saveDir);
    end
    save([saveDir,'AngDevVars.mat'],'AD')
else
    load([saveDir,'AngDevVars.mat'])
end


% Plot weigthed mean angular deviation variables per condition
% --------------------------- User Input ----------------------------------
% Here you can specify which angular deviation related data you'd like to
% plot for visualization. All data weighted by total bout length per fly.
% 'Mean' will plot the mean angular deviation per fly
% 'Abs' will plot the mean absolute angular deviation per fly
% 'Binary' will plot the binarized angular deviation per fly
% 'Std' will plot the standard deviation of mean angular deviation per fly
% 'Bout' will plot the mean forward bout lengths per fly
% Accepted arguments: ['Mean','Abs','Binary','Std','Bout']
AngDevToPlot = 'Mean';

% Here you can filter the data plotted based on the mean forward velocity
% of the fly during a forward bout. As of 221130, definitions of fast and
% slow bouts are static and pre-determined during data loading phase.
% 'Slow' and 'Fast' will plot only the slow and fast bouts, 'All' will plot
% all the data. Note that Slow + Fast does not have to equal All bouts.
% Accepted arguments: ['Slow','Fast','All']
SpeedToPlot = 'All';

sig = 1; % 1 if you want to overlay significance stars (MWW non-param test)
fs = 14; % font size for plot labels
alpha = 0.5; % transparency of the dots
ms = 32; % dot size of the mean data point
scaleFactor = 0.1; % scale the dot size
scatterFactor = 0.26; % add scatter (fraction between 0 and 1) to the dots
limsX = [0 size(AD.Wmu,2)+1]; % axis limits

%----------------------------end of user input-----------------------------
%AD.seq = seq'; % add stimulus info (no 'Dark' for Light-only analysis)
%AD.seq = horzcat({'Dark'},AD.seq'); % add stimulus info


AD.seq = [{'Dark'}, AD.seq(:)'];
switch size(AD.seq,2) % change plotting order so 10deg and RG comes last
    case 1 % if there is only dark
        conds = 1; xTicks = 1;
    case 2
        % 1 + 10 NG  (bIPS Translational)
        conds = [2,1]; xTicks = 1:size(AD.Wmu,2);
    case 3
        % Dark + 1 + 10 NG (DNp15 Ricin A)
        conds = [1,3,2]; xTicks = 1:size(AD.Wmu,2);
    case 7        
        % Dark + 1 + 5 + 10 NG|RG  (DNp15 hsFLP)
        conds = [1,4,6,2,5,7,3]; xTicks = 1:size(AD.Wmu,2);
    case 8
        conds = [3,5,7,1,4,6,8,2]; xTicks = 1:8;
    case 9 % assumes dark is present
        % Dark + 1 + 5 + 10 NG|RG  (bIPS hsFLP)
        % conds = [1,4,6,8,2,5,7,9,3]; xTicks = 1:size(AD.Wmu,2);
        conds = [1,4,2]; xTicks = 1:3;     
end

% calculate maximum y values across data points to have the significance
% star plots well separated over the data points
my = zeros(1,length(conds));

f=figure('WindowState','maximized');

i = find(contains(AD.seq,'10NG'));
sf = scaleFactor;
hold on;

if contains(behDir,'DNp15')
    col = cReg(1,:);
elseif contains(behDir,'bIPS')
    col = cReg(2,:);
else
    col = [0,0,0];
end

for g = 1:length(genotypes)

    switch AngDevToPlot
        % vp: variable to plot | yl: y-axis label | sn: save name |
        case 'Mean'
            switch SpeedToPlot
                case 'Slow'
                    vp = AD.WmuS;
                case 'Fast'
                    vp = AD.WmuF;
                case 'All'
                    vp = AD.Wmu;
                otherwise
                    error('Unknown/unspecified SpeedToPlot variable!')
            end
            sn = 'AngDevPerCond';
            yl = 'Ang. Dev. (deg/mm)';
            limsY = [-7 7];
        case 'Abs'
            switch SpeedToPlot
                case 'Slow'
                    vp = AD.AbsWmuS;
                case 'Fast'
                    vp = AD.AbsWmuF;
                case 'All'
                    vp = AD.AbsWmu;
                otherwise
                    error('Unknown/unspecified SpeedToPlot variable!')
            end
            sn = 'AbsAngDevPerCond';
            yl = 'Abs. Ang. Dev. (deg/mm)';
            limsY = [0 5];
        case 'Std'
            switch SpeedToPlot
                case 'Slow'
                    vp = AD.WstdS;
                case 'Fast'
                    vp = AD.WstdF;
                case 'All'
                    vp = AD.Wstd;
                otherwise
                    error('Unknown/unspecified SpeedToPlot variable!')
            end
            sn = 'StdAngDevPerCond';
            yl = 'std of Ang. Dev. (a.u.)';
            limsY = [0 10];
        case 'Bout'
            switch SpeedToPlot
                case 'Slow'
                    vp = AD.BlmuS;
                case 'Fast'
                    vp = AD.BlmuF;
                case 'All'
                    vp = AD.Blmu;
                otherwise
                    error('Unknown/unspecified SpeedToPlot variable!')
            end
            sn = 'FwBoutLengthPerCond';
            yl = 'Fw bout length (frames)';
            limsY = [20 80];
        case 'Binary'
            switch SpeedToPlot
                case 'Slow'
                    error('Binarized plots for Slow not implemented!')
                case 'Fast'
                    error('Binarized plots for Fast not implemented!')
                case 'All'
                    vp = AD.BWmu;
                otherwise
                    error('Unknown/unspecified SpeedToPlot variable!')
            end
            sn = 'AngBiasPerCond';
            yl = 'Angular Bias';
            limsY = [-1.1 1.1];
        otherwise
            error(['Don''t know what to plot, ',...
                'please check the AngDevToPlot variable!'])
    end
    switch SpeedToPlot
        % bl: bout length data used to weigh the per fly data
        case 'Slow'
            bl = AD.BlSumS;
        case 'Fast'
            bl = AD.BlSumF;
        case 'All'
            bl = AD.BlSum;
        otherwise
            error('Unknown/unspecified SpeedToPlot variable!')
    end

    % calculate weighted grand mean and its variance & std deviation
    Wmean = nansum((vp{1,gi(g)}(:,i).*...
        bl{1,gi(g)}(:,i))/nansum(bl{1,gi(g)}(:,i)));
    Wvar = nansum((((vp{1,gi(g)}(:,i) - Wmean).^2).*...
        bl{1,gi(g)} (:,i))/nansum(bl{1,gi(g)}(:,i)));
    Wstd = sqrt(Wvar);
    % define x-axis location of each dot
    sc = (rand(size(AD.Blmu{1,gi(g)},1),1)-0.5)*scatterFactor;
    x = (ones(size(AD.Blmu{1,gi(g)},1),1)*g) + sc;
    % plot the data, add mean and std error bars
    scatter(x,vp{1,gi(g)}(:,i),...
        (bl{1,gi(g)}(:,i)*sf),...
        'MarkerFaceColor',col,'MarkerEdgeColor','none',...
        'MarkerFaceAlpha',alpha);
    plot(g,Wmean,'k.','MarkerSize',ms)
    nf = size(AD.Blmu{1,gi(g)},1); % number of flies in this genotype
    eh = errorbar(g,Wmean,Wstd/sqrt(nf));
    set(eh,'LineWidth',2,'Color','k')
    maxbl = nanmax(vp{1,gi(g)}(:,i));
    my(i) = max([my(i),maxbl,limsY(2)*0.4]);
    my(i) = min(my(i),limsY(2)*0.4);
end
if sig
    % Perform Mann-Whitney-Wilcoxon non parametric test and plot
    % significance indicators if p<0.05
    nck = nchoosek(1:g,2);
    %         % bIPS TrpA1 Kir combined plots only: remove unnecessary comparisons
    %         nck = nck([1,10:end],:);
    nph = size(nck,1); % number of comparisons made (used for adjusting the
    % minimum P value required for significance - a.k.a. Bonferroni correction)
    for ii = 1 : size(nck,1)
        pair1 = vp{1,gi(nck(ii,1))}(:,i);
        pair2 = vp{1,gi(nck(ii,2))}(:,i);
        pair1(isnan(pair1)) = []; pair2(isnan(pair2)) = [];
        stats{ii} = mwwtest(pair1,pair2,0);
        % uncomment below if you want to display all p-values per pair
        % disp([AD.seq{i},'-',gTypes{nck(ii,1)},'-',gTypes{nck(ii,2)},'-',num2str(stats{ii}.p(2))])
        if stats{ii}.p < 0.05/nph
            % uncomment below if you want to display significant p-values per pair
            %                 disp([AD.seq{i},'-',gTypes{nck(ii,1)},'-',gTypes{nck(ii,2)},'-',num2str(stats{ii}.p(2))])
            sigline([(nck(ii,1)),(nck(ii,2))],num2str(stats{ii}.p(2)),my(i)+my(i)*0.1*ii);
        end
    end
end
% formatting the plot
% set axis limits
xlim(limsX); ylim(limsY);
% add tick and axis labels
ylabel(yl)
set(gca,'XTick',xTicks,'XTickLabel',gOrder,'FontSize',fs,...
    'XTickLabelRotation',45)
if contains(AngDevToPlot,'Mean') || contains(AngDevToPlot,'Binary')
    % add dashed line at 0
    line(limsX,[0 0],'Color','k','LineWidth',0.5,'LineStyle','--')
end
% text overlay of speed
xL = 2; % text location on the x-axis
yL = my(i)-my(i)*0.11*ii; % text location on the y-axis
switch SpeedToPlot
    case 'Slow'
        T = text(xL,yL,['Slow (<',num2str(AD.VfwThr(1)),'mm/s)']);
    case 'Fast'
        T = text(xL,yL,['Fast (>',num2str(AD.VfwThr(2)),'mm/s)']);
end


%% Figure 7e

% bIPS hsFlp Kir subfolder
behDir = 'bIPSKir\';
cond = {'Dark','Light'}; % experimental conditions
genotypes = {'None','Left','Right','Both'};
gTypes = {'None','Left','Right','Both'};
% save directory
saveDir = [baseDir,'Analysis\bIPSKir\'];
[~,gi] = ismember(gOrder,gTypes); % plotting order indices

% load angular deviation data. If AD.mat doesn't exist, calculate angular
% deviations, store in AD.mat and save it on disk
if ~exist([saveDir,'AngDevVars.mat'],'file')
    % initialize the structures. Each structure will have as many columns
    % as there are genotypes. Cell arrays within each coulm will have size
    % of MxN where M is the number of flies and N is number of total
    % conditions (dark + 8 different light conditions as of 220131).
    % Sequence of light conditions will be stored in the seq variable below
    lightCondNum = 2;
    VfwThrS = 14; % mm/s | the maximum forward velocity that will be
    % considered as slow walking
    VfwThrF = 18; % mm/s | the minimum forward velocity that will be
    % considered as fast walking
    AD = struct();
    AD.VfwThr = [VfwThrS,VfwThrF];
    % structure fields for all forward bouts
    AD.Wmu = cell(1,length(genotypes)); % weighted mean
    AD.Wstd = cell(1,length(genotypes)); % std of the weighted mean
    AD.AbsWmu = cell(1,length(genotypes)); % weighted mean of abs angdev
    AD.BWmu = cell(1,length(genotypes)); %weighted mean of binarized angdev
    AD.BlSum = cell(1,length(genotypes)); % bout length summed per fly
    AD.Bl = cell(1,length(genotypes)); % bout length per bout per fly
    AD.Adev = cell(1,length(genotypes)); % bout length per bout per fly
    AD.Blmu = cell(1,length(genotypes)); % mean bout length per fly
    AD.Vfmu = cell(1,length(genotypes)); % mean Vfw per bout per fly
    AD.VFperGen = cell(1,length(genotypes)); % Vfw values per genotype
    % structure fields for "fast" forward bouts
    AD.WmuF = cell(1,length(genotypes)); % weighted mean
    AD.WstdF = cell(1,length(genotypes)); % std of the weighted mean
    AD.AbsWmuF = cell(1,length(genotypes)); % weighted mean of abs angdev
    AD.BlSumF = cell(1,length(genotypes)); % bout length summed per fly
    AD.BlF = cell(1,length(genotypes)); % bout length per bout per fly
    AD.BlmuF = cell(1,length(genotypes)); % mean bout length per fly
    AD.VfmuF = cell(1,length(genotypes)); % mean Vfw per bout per fly
    % structure fields for "slow" forward bouts
    AD.WmuS = cell(1,length(genotypes)); % weighted mean
    AD.WstdS = cell(1,length(genotypes)); % std of the weighted mean
    AD.AbsWmuS = cell(1,length(genotypes)); % weighted mean of abs angdev
    AD.BlSumS = cell(1,length(genotypes)); % bout length summed per fly
    AD.BlS = cell(1,length(genotypes)); % bout length per bout per fly
    AD.BlmuS = cell(1,length(genotypes)); % mean bout length per fly
    AD.VfmuS = cell(1,length(genotypes)); % mean Vfw per bout per fly
    % All angular deviation values concatenated per genotype
    AD.AdevAll = cell(1,length(genotypes));

    for c = 1:length(cond)
        for g = 1:length(genotypes) % for each genotype
        
        % directory where the files are located
        path = [baseDir,behDir,cond{c},'\',genotypes{g},'\'];
        disp([behDir,cond{c},'\',genotypes{g},'\'])
        % load the low resolution data relative to the forward segments
        % together with the stimuli sequence
        [VForwSeg,seq] = GetForwSegVData(path, params);
        % NOTE: Fw speed cutoff is stored in params.vft and it can be
        % changed based on slow/fast forward speed analysis preference

        % binning windows for the angular displacements (deg/mm)
        DevCents = [-1000 -3.75 -2.5 -1.25 0 1.25 2.5 3.75 1000];

        % initialize variables
        VR = cell(length(DevCents)-1,1); % rotational velocity
        VF = cell(length(DevCents)-1,1); % forward velocity
        VS = cell(length(DevCents)-1,1); % side velocity
        NB = cell(length(DevCents)-1,1); % number of bouts
        angDevAll = [];
        % Length of each forward walking bout per fly
        Blength = cell(size(VForwSeg,2),size(VForwSeg,1));
        % Angular deviation per bout per fly
        Adev = cell(size(VForwSeg,2),size(VForwSeg,1)); % signed
        AdevAbs = cell(size(VForwSeg,2),size(VForwSeg,1)); % absolute
        AdevB = cell(size(VForwSeg,2),size(VForwSeg,1)); % binarized
        % Mean forward velocity per bout per fly
        Vfw = cell(size(VForwSeg,2),size(VForwSeg,1));

        for jk = 1:length(DevCents)-1 % for each angular deviation bin
            % define boundaries of each bin
            minAngDev = DevCents(jk);
            maxAngDev = DevCents(jk+1);

            for s = 1:size(VForwSeg,1) % for each stimulus
                for n = 1:size(VForwSeg,2) % for each fly
                    for i = 1:length(VForwSeg{s,n}) % for each fw bout
                        if ~isempty(VForwSeg{s,n}{i})% if there is data
                            % low res speed data for this fw segment
                            vr = VForwSeg{s,n}{i}.VrLR; %Vrotational
                            vf = VForwSeg{s,n}{i}.VfLR; %Vforward
                            vs = VForwSeg{s,n}{i}.VsLR; %Vsideways
                            vt = sqrt(vf.*vf + vs.*vs); %Vtranslational
                            angDev = mean(vr./vt); % angular deviation
                            angDevB = angDev/abs(angDev); % binarized

                            %if the angular deviation belong to the bin
                            if angDev > minAngDev && angDev < maxAngDev
                                % store the data appropriately
                                VR{jk} = vertcat(VR{jk}, vr);
                                VF{jk} = vertcat(VF{jk}, vf);
                                VS{jk} = vertcat(VS{jk}, vs);
                                NB{jk} = vertcat(NB{jk}, length(vs));
                                angDevAll = vertcat(angDevAll,angDev);
                                Adev{n,s} = vertcat(Adev{n,s},angDev);
                                AdevB{n,s} = vertcat(AdevB{n,s},angDevB);
                                AdevAbs{n,s} = vertcat(AdevAbs{n,s},...
                                    abs(angDev));
                                Vfw{n,s} = vertcat(Vfw{n,s},mean(vf));
                                Blength{n,s} = vertcat(Blength{n,s},...
                                    length(vs));
                            end
                        end
                    end
                end
            end
        end

        % initialize variables
        % Angular deviation per fly weighted by ALL fw bout length
        AngDevWmu = nan(size(Adev));
        AngDevWvar = nan(size(Adev));
        AngDevAbsWmu = nan(size(Adev));
        AngDevBWmu = nan(size(AdevB));
        BlSum = nan(size(Adev));
        BlW = nan(size(Adev));
        AD.AdevAll{g} = vertcat(AD.AdevAll{g},angDevAll);
        % Angular deviation per fly weighted by FAST fw bout length
        AngDevWmuF = nan(size(Adev));
        AngDevWvarF = nan(size(Adev));
        AngDevAbsWmuF = nan(size(Adev));
        BlSumF = nan(size(Adev));
        AdevF = cell(size(Adev));
        AdevAbsF = cell(size(AdevAbs));
        BlengthF = cell(size(Blength));
        % Angular deviation per fly weighted by SLOW fw bout length
        AngDevWmuS = nan(size(Adev));
        AngDevWvarS = nan(size(Adev));
        AngDevAbsWmuS = nan(size(Adev));
        BlSumS = nan(size(Adev));
        AdevS = cell(size(Adev));
        BlengthS = cell(size(Blength));
        AdevAbsS = cell(size(AdevAbs));

        % logical index arrays for slow and fast bouts
        iFast =cellfun(@(x) x >= VfwThrF, Vfw, 'UniformOutput', false);
        iSlow =cellfun(@(x) x < VfwThrS, Vfw, 'UniformOutput', false);

        for n = 1:size(Adev,1) % for each fly
            for i = 1:size(Adev,2) % for each stimulus

                % Angular deviations of all bouts calculated as a
                % weighted mean based on bout lengths
                AngDevWmu(n,i)=...
                    nansum((Adev{n,i}.*Blength{n,i})/...
                    nansum(Blength{n,i}));
                AngDevAbsWmu(n,i)=...
                    nansum((AdevAbs{n,i}.*Blength{n,i})/...
                    nansum(Blength{n,i}));
                AngDevBWmu(n,i)=...
                    nansum((AdevB{n,i}.*Blength{n,i})/...
                    nansum(Blength{n,i}));
                AngDevWvar(n,i)=...
                    nansum((((Adev{n,i}-AngDevWmu(n,i)).^2)...
                    .*Blength{n,i})/nansum(Blength{n,i}));
                BlSum(n,i) = nansum(Blength{n,i});

                % Angular deviations of fast bouts calculated as a
                % weighted mean based on bout lengths
                AdevF{n,i} = Adev{n,i}(iFast{n,i}); % fast bouts only
                AdevAbsF{n,i} = AdevAbs{n,i}(iFast{n,i});
                BlengthF{n,i} = Blength{n,i}(iFast{n,i});
                AngDevWmuF(n,i)=...
                    nansum((AdevF{n,i}.*BlengthF{n,i})/...
                    nansum(BlengthF{n,i}));
                AngDevAbsWmuF(n,i)=...
                    nansum((AdevAbsF{n,i}.*BlengthF{n,i})/...
                    nansum(BlengthF{n,i}));
                AngDevWvarF(n,i)=...
                    nansum((((AdevF{n,i}-AngDevWmuF(n,i)).^2)...
                    .*BlengthF{n,i})/nansum(BlengthF{n,i}));
                BlSumF(n,i) = nansum(BlengthF{n,i});

                % Angular deviations of slow bouts calculated as a
                % weighted mean based on bout lengths
                AdevS{n,i} = Adev{n,i}(iSlow{n,i}); % slow bouts only
                AdevAbsS{n,i} = AdevAbs{n,i}(iSlow{n,i});
                BlengthS{n,i} = Blength{n,i}(iSlow{n,i});
                AngDevWmuS(n,i)=...
                    nansum((AdevS{n,i}.*BlengthS{n,i})/...
                    nansum(BlengthS{n,i}));
                AngDevAbsWmuS(n,i)=...
                    nansum((AdevAbsS{n,i}.*BlengthS{n,i})/...
                    nansum(BlengthS{n,i}));
                AngDevWvarS(n,i)=...
                    nansum((((AdevS{n,i}-AngDevWmuS(n,i)).^2)...
                    .*BlengthS{n,i})/nansum(BlengthS{n,i}));
                BlSumS(n,i) = nansum(BlengthS{n,i});
            end
        end

        % remove zeros from the data
        AngDevWmu(BlSum==0) = nan;
        AngDevWvar(BlSum==0) = nan;
        AngDevAbsWmu(BlSum==0) = nan;
        AngDevBWmu(BlSum==0) = nan;
        BlSum(BlSum==0) = nan;
        AngDevWmuF(BlSumF==0) = nan;
        AngDevWvarF(BlSumF==0) = nan;
        AngDevAbsWmuF(BlSumF==0) = nan;
        BlSumF(BlSumF==0) = nan;
        AngDevWmuS(BlSumS==0) = nan;
        AngDevWvarS(BlSumS==0) = nan;
        AngDevAbsWmuS(BlSumS==0) = nan;
        BlSumS(BlSumS==0) = nan;
        % standard deviation of the weighted means
        AngDevWstd = sqrt(AngDevWvar);
        AngDevWstdF = sqrt(AngDevWvarF);
        AngDevWstdS = sqrt(AngDevWvarS);

        % horzcat concatenates in the opposite of what we want when the
        % initial array is empty, so we swap the positions of variables
        % during the first concatenation (which happens in 'Dark'
        % condition). If you get errors here, make sure that the first
        % element of conds variable is 'Dark'.
        switch cond{c}
            case 'Dark'
                AD.Wmu{g} = horzcat(AngDevWmu,AD.Wmu{g});
                AD.Wstd{g} = horzcat(AngDevWstd,AD.Wstd{g});
                AD.AbsWmu{g} = horzcat(AngDevAbsWmu,AD.AbsWmu{g});
                AD.BWmu{g} = horzcat(AngDevBWmu,AD.BWmu{g});
                AD.BlSum{g} = horzcat(BlSum,AD.BlSum{g});
                AD.Bl{g} = horzcat(Blength,AD.Bl{g});
                AD.Adev{g} = horzcat(Adev,AD.Adev{g});
                AD.WmuF{g} = horzcat(AngDevWmuF,AD.WmuF{g});
                AD.WstdF{g} = horzcat(AngDevWstdF,AD.WstdF{g});
                AD.AbsWmuF{g} = horzcat(AngDevAbsWmuF,AD.AbsWmuF{g});
                AD.BlSumF{g} = horzcat(BlSumF,AD.BlSumF{g});
                AD.BlF{g} = horzcat(BlengthF,AD.BlF{g});
                AD.WmuS{g} = horzcat(AngDevWmuS,AD.WmuS{g});
                AD.WstdS{g} = horzcat(AngDevWstdS,AD.WstdS{g});
                AD.AbsWmuS{g} = horzcat(AngDevAbsWmuS,AD.AbsWmuS{g});
                AD.BlSumS{g} = horzcat(BlSumS,AD.BlSumS{g});
                AD.BlS{g} = horzcat(BlengthS,AD.BlS{g});

            case {'Light','Light Bilateral'}
                AD.Wmu{g} = horzcat(AD.Wmu{g},AngDevWmu);
                AD.Wstd{g} = horzcat(AD.Wstd{g},AngDevWstd);
                AD.AbsWmu{g} = horzcat(AD.AbsWmu{g},AngDevAbsWmu);
                AD.BWmu{g} = horzcat(AD.BWmu{g},AngDevBWmu);
                AD.BlSum{g} = horzcat(AD.BlSum{g},BlSum);
                AD.Bl{g} = horzcat(AD.Bl{g},Blength);
                AD.Adev{g} = horzcat(AD.Adev{g},Adev);
                AD.WmuF{g} = horzcat(AD.WmuF{g},AngDevWmuF);
                AD.WstdF{g} = horzcat(AD.WstdF{g},AngDevWstdF);
                AD.AbsWmuF{g} = horzcat(AD.AbsWmuF{g},AngDevAbsWmuF);
                AD.BlSumF{g} = horzcat(AD.BlSumF{g},BlSumF);
                AD.BlF{g} = horzcat(AD.BlF{g},BlengthF);
                AD.WmuS{g} = horzcat(AD.WmuS{g},AngDevWmuS);
                AD.WstdS{g} = horzcat(AD.WstdS{g},AngDevWstdS);
                AD.AbsWmuS{g} = horzcat(AD.AbsWmuS{g},AngDevAbsWmuS);
                AD.BlSumS{g} = horzcat(AD.BlSumS{g},BlSumS);
                AD.BlS{g} = horzcat(AD.BlS{g},BlengthS);
                AD.seq = horzcat({'Dark'},seq'); % add stimulus info
        end

        %% Figure 7b - Plot example trajectories
        xT = cell(1,1);
        yT = cell(1,1);
        indd = cell(1,1);
        a = 1;
        cmap1 = jet(length(VR));
        xTall = [];
        yTall = [];
        % iterate through the angular deviation bins
        for i = 2 : length(NB)-1
            nbs = NB{i};
            nbScs = [1; cumsum(nbs)];
            % iterate through the forward segments
            for j = 1 : length(nbs)
                vr = VR{i}(nbScs(j):nbScs(j+1))/60;
                vf = VF{i}(nbScs(j):nbScs(j+1))/60;
                vs = VS{i}(nbScs(j):nbScs(j+1))/60;
                x = [];
                y = [];
                th = [];
                if length(vr) < 180
                    % get centered X and Y
                    for ij = 1 : length(vr)
                        if isempty(x)
                            th = vertcat(th,pi*vr(ij)/180);
                            x = vertcat(x,vf(ij)*sin(0)+vs(ij)*cos(0));
                            y = vertcat(y,vf(ij)*cos(0)-vs(ij)*sin(0));
                        else
                            th = vertcat(th,th(end)+pi*vr(ij)/180);
                            x = vertcat(x,x(end)+...
                                vf(ij)*sin(th(ij-1))...
                                +vs(ij)*cos(th(ij-1)));
                            y = vertcat(y,y(end)+...
                                vf(ij)*cos(th(ij-1))...
                                -vs(ij)*sin(th(ij-1)));
                        end
                    end
                end
                xT{a} = x;
                yT{a} = y;
                indd{a} = i;
                if a == 1
                    xTall = x';
                    yTall = y';
                else
                    if length(x) > size(xTall,2)
                        xTall(:,(size(xTall,2)+1):length(x)) = NaN;
                        yTall(:,(size(yTall,2)+1):length(y)) = NaN;
                        xTall(a,:) = x';
                        yTall(a,:) = y';
                    else
                        xTall(a,:) = NaN;
                        yTall(a,:) = NaN;
                        xTall(a,1:length(x)) = x';
                        yTall(a,1:length(y)) = y';
                    end
                end
                a = a + 1;
            end
        end
        yTallBin = 0:round(max(nanmax(yTall)));
        xTallBin = nan(1,length(yTallBin));
        yTallR = round(yTall);
        for i = 1:(round(max(nanmax(yTall)))+1)
            xTallBin(i)=nanmean(xTall(yTallR==i-1));
        end
        % plot a color coded sample of the forward segments
        inds = randperm(length(xT));
        frac = 100; % percentage of data to be plotted
        ythr = 10; % threshold of forward translation (mm in y axis)
        inds = inds(1:round(length(xT)*frac/100));

        fs = 18; % font size
        f=figure('WindowState','maximized');
        hold on
        for i = 1 : length(inds)
            if max(yT{inds(i)})>ythr
                plot(yT{inds(i)}, xT{inds(i)},...
                    'color', cmap1(indd{inds(i)},:),'LineWidth',1.5)
            end
        end
        plot(yTallBin,xTallBin,'color','k','LineWidth',2)

        axis([-2 30 -25 25])
        limsX = xlim;
        line(limsX,[0 0],'Color','k','LineWidth',0.5,'LineStyle','--')
        %             title(['bIPS-',gTypes{g}]) % add title, tick and axis labels
        %             title(['DNp15-',genotypes{g}]) % add title, tick and axis labels
        title(['DNp15-RicinA-',gTypes{g}]) % add title, tick and axis labels
        ylabel('X position (mm)'); xlabel('Y position (mm)')
        set(gca,'FontSize',fs)
        caxis([DevCents(2) DevCents(end-1)])
        colormap(cmap1)
        ch = colorbar('Ticks',DevCents(2:end-1));
        ch.Label.String = 'Ang. Dev. (deg/mm)';
        ch.Label.Rotation = 270;
        ts = timestamp(); % get current date
        if SaveFigures
            set(gcf,'renderer','Painters')
            print(f,strcat(saveDir,ts,'-Paths-',...
                genotypes{g},'-',cond{c},'-',num2str(frac),...
                'perc-yThr',num2str(ythr)),'-dtiff') %#ok<*UNRCH>
            print(f,strcat(saveDir,ts,'-Paths-',...
                genotypes{g},'-',cond{c},'-',num2str(frac),...
                'perc-yThr',num2str(ythr)),'-depsc')
            close(f)
        end
        end
    end

    %     AD.AdevAll = AdevAll;
    for g = 1:length(genotypes) % for each genotype
        nFly = size(AD.Bl{1,g}(:,1),1);
        Blmean = nan(nFly,size(AD.seq,2));
        BlmeanF = nan(nFly,size(AD.seq,2));
        BlmeanS = nan(nFly,size(AD.seq,2));
        for i = 1:size(AD.seq,2) % for each condition
            for n = 1:nFly % for each fly
                Blmean(n,i) = nanmean(AD.Bl{1,g}{n,i});
                BlmeanF(n,i) = nanmean(AD.BlF{1,g}{n,i});
                BlmeanS(n,i) = nanmean(AD.BlS{1,g}{n,i});
            end
        end
        AD.Blmu{1,g} = Blmean;
        AD.BlmuF{1,g} = BlmeanF;
        AD.BlmuS{1,g} = BlmeanS;
    end
    if ~exist(saveDir, 'dir')
        mkdir(saveDir);
    end
    save([saveDir,'AngDevVars.mat'],'AD')
else
    load([saveDir,'AngDevVars.mat'])
end


% Plot weigthed mean angular deviation variables per condition
% --------------------------- User Input ----------------------------------
% Here you can specify which angular deviation related data you'd like to
% plot for visualization. All data weighted by total bout length per fly.
% 'Mean' will plot the mean angular deviation per fly
% 'Abs' will plot the mean absolute angular deviation per fly
% 'Binary' will plot the binarized angular deviation per fly
% 'Std' will plot the standard deviation of mean angular deviation per fly
% 'Bout' will plot the mean forward bout lengths per fly
% Accepted arguments: ['Mean','Abs','Binary','Std','Bout']
AngDevToPlot = 'Mean';

% Here you can filter the data plotted based on the mean forward velocity
% of the fly during a forward bout. As of 221130, definitions of fast and
% slow bouts are static and pre-determined during data loading phase.
% 'Slow' and 'Fast' will plot only the slow and fast bouts, 'All' will plot
% all the data. Note that Slow + Fast does not have to equal All bouts.
% Accepted arguments: ['Slow','Fast','All']
SpeedToPlot = 'All';

sig = 1; % 1 if you want to overlay significance stars (MWW non-param test)
fs = 14; % font size for plot labels
alpha = 0.5; % transparency of the dots
ms = 32; % dot size of the mean data point
scaleFactor = 0.1; % scale the dot size
scatterFactor = 0.26; % add scatter (fraction between 0 and 1) to the dots
limsX = [0 size(AD.Wmu,2)+1]; % axis limits

%----------------------------end of user input-----------------------------

switch size(AD.seq,2) % change plotting order so 10deg and RG comes last
    case 1 % if there is only dark
        conds = 1; xTicks = 1; colors = [202,202,202]/256;
    case 2
        % 1 + 10 NG  (bIPS Translational)
        conds = [2,1]; xTicks = 1:size(AD.Wmu,2);
    case 3
        % Dark + 1 + 10 NG (DNp15 Ricin A)
        conds = [1,3,2]; xTicks = 1:size(AD.Wmu,2);
    case 7        
        % Dark + 1 + 5 + 10 NG|RG  (DNp15 hsFLP)
        conds = [1,4,6,2,5,7,3]; xTicks = 1:size(AD.Wmu,2);
    case 8
        conds = [3,5,7,1,4,6,8,2]; xTicks = 1:8;
    case 9 % assumes dark is present
        % Dark + 1 + 5 + 10 NG|RG  (bIPS hsFLP)
        conds = [1,4,6,8,2,5,7,9,3]; xTicks = 1:size(AD.Wmu,2);
        % conds = [1,4,2]; xTicks = 1:3; ci=[1,2,4];
end

% calculate maximum y values across data points to have the significance
% star plots well separated over the data points
my = zeros(1,length(conds));

f=figure('WindowState','maximized');

i = find(contains(AD.seq,'10NG'));
sf = scaleFactor;
hold on;

if contains(behDir,'DNp15')
    col = cReg(1,:);
elseif contains(behDir,'bIPS')
    col = cReg(2,:);
else
    col = [0,0,0];
end

for g = 1:length(genotypes)

    switch AngDevToPlot
        % vp: variable to plot | yl: y-axis label | sn: save name |
        case 'Mean'
            switch SpeedToPlot
                case 'Slow'
                    vp = AD.WmuS;
                case 'Fast'
                    vp = AD.WmuF;
                case 'All'
                    vp = AD.Wmu;
                otherwise
                    error('Unknown/unspecified SpeedToPlot variable!')
            end
            sn = 'AngDevPerCond';
            yl = 'Ang. Dev. (deg/mm)';
            limsY = [-7 7];
        case 'Abs'
            switch SpeedToPlot
                case 'Slow'
                    vp = AD.AbsWmuS;
                case 'Fast'
                    vp = AD.AbsWmuF;
                case 'All'
                    vp = AD.AbsWmu;
                otherwise
                    error('Unknown/unspecified SpeedToPlot variable!')
            end
            sn = 'AbsAngDevPerCond';
            yl = 'Abs. Ang. Dev. (deg/mm)';
            limsY = [0 5];
        case 'Std'
            switch SpeedToPlot
                case 'Slow'
                    vp = AD.WstdS;
                case 'Fast'
                    vp = AD.WstdF;
                case 'All'
                    vp = AD.Wstd;
                otherwise
                    error('Unknown/unspecified SpeedToPlot variable!')
            end
            sn = 'StdAngDevPerCond';
            yl = 'std of Ang. Dev. (a.u.)';
            limsY = [0 10];
        case 'Bout'
            switch SpeedToPlot
                case 'Slow'
                    vp = AD.BlmuS;
                case 'Fast'
                    vp = AD.BlmuF;
                case 'All'
                    vp = AD.Blmu;
                otherwise
                    error('Unknown/unspecified SpeedToPlot variable!')
            end
            sn = 'FwBoutLengthPerCond';
            yl = 'Fw bout length (frames)';
            limsY = [20 80];
        case 'Binary'
            switch SpeedToPlot
                case 'Slow'
                    error('Binarized plots for Slow not implemented!')
                case 'Fast'
                    error('Binarized plots for Fast not implemented!')
                case 'All'
                    vp = AD.BWmu;
                otherwise
                    error('Unknown/unspecified SpeedToPlot variable!')
            end
            sn = 'AngBiasPerCond';
            yl = 'Angular Bias';
            limsY = [-1.1 1.1];
        otherwise
            error(['Don''t know what to plot, ',...
                'please check the AngDevToPlot variable!'])
    end
    switch SpeedToPlot
        % bl: bout length data used to weigh the per fly data
        case 'Slow'
            bl = AD.BlSumS;
        case 'Fast'
            bl = AD.BlSumF;
        case 'All'
            bl = AD.BlSum;
        otherwise
            error('Unknown/unspecified SpeedToPlot variable!')
    end

    % calculate weighted grand mean and its variance & std deviation
    Wmean = nansum((vp{1,gi(g)}(:,i).*...
        bl{1,gi(g)}(:,i))/nansum(bl{1,gi(g)}(:,i)));
    Wvar = nansum((((vp{1,gi(g)}(:,i) - Wmean).^2).*...
        bl{1,gi(g)} (:,i))/nansum(bl{1,gi(g)}(:,i)));
    Wstd = sqrt(Wvar);
    % define x-axis location of each dot
    sc = (rand(size(AD.Blmu{1,gi(g)},1),1)-0.5)*scatterFactor;
    x = (ones(size(AD.Blmu{1,gi(g)},1),1)*g) + sc;
    % plot the data, add mean and std error bars
    scatter(x,vp{1,gi(g)}(:,i),...
        (bl{1,gi(g)}(:,i)*sf),...
        'MarkerFaceColor',col,'MarkerEdgeColor','none',...
        'MarkerFaceAlpha',alpha);
    plot(g,Wmean,'k.','MarkerSize',ms)
    nf = size(AD.Blmu{1,gi(g)},1); % number of flies in this genotype
    eh = errorbar(g,Wmean,Wstd/sqrt(nf));
    set(eh,'LineWidth',2,'Color','k')
    maxbl = nanmax(vp{1,gi(g)}(:,i));
    my(i) = max([my(i),maxbl,limsY(2)*0.4]);
    my(i) = min(my(i),limsY(2)*0.4);
end
if sig
    % Perform Mann-Whitney-Wilcoxon non parametric test and plot
    % significance indicators if p<0.05
    nck = nchoosek(1:g,2);
    nph = size(nck,1); % number of comparisons made (used for adjusting the
    % minimum P value required for significance - a.k.a. Bonferroni correction)
    for ii = 1 : size(nck,1)
        pair1 = vp{1,gi(nck(ii,1))}(:,i);
        pair2 = vp{1,gi(nck(ii,2))}(:,i);
        pair1(isnan(pair1)) = []; pair2(isnan(pair2)) = [];
        stats{ii} = mwwtest(pair1,pair2,0);
        if stats{ii}.p < 0.05/nph
           sigline([(nck(ii,1)),(nck(ii,2))],num2str(stats{ii}.p(2)),my(i)+my(i)*0.1*ii);
        end
    end
end
% formatting the plot
% set axis limits
xlim(limsX); ylim(limsY);
% add tick and axis labels
ylabel(yl)
set(gca,'XTick',xTicks,'XTickLabel',gOrder,'FontSize',fs,...
    'XTickLabelRotation',45)
if contains(AngDevToPlot,'Mean') || contains(AngDevToPlot,'Binary')
    % add dashed line at 0
    line(limsX,[0 0],'Color','k','LineWidth',0.5,'LineStyle','--')
end
% text overlay of speed
xL = 2; % text location on the x-axis
yL = my(i)-my(i)*0.11*ii; % text location on the y-axis
switch SpeedToPlot
    case 'Slow'
        T = text(xL,yL,['Slow (<',num2str(AD.VfwThr(1)),'mm/s)']);
    case 'Fast'
        T = text(xL,yL,['Fast (>',num2str(AD.VfwThr(2)),'mm/s)']);
end

%% Extended Data Figure 10a and 10b

% % DNp15 hsFlp RicinA
behDir = 'DNp15Ricin\';

cond = {'Dark','Light'}; % experimental conditions
genotypes = {'None','Left','Right','Both'};
gTypes = {'None','Left','Right','Both'};

saveDir = [baseDir,'Analysis\DNp15Ricin\']; % save directory

[~,gi] = ismember(gOrder,gTypes); % plotting order indices


% load angular deviation data. If AD.mat doesn't exist, calculate angular
% deviations, store in AD.mat and save it on disk
if ~exist([saveDir,'AngDevVars.mat'],'file')
    % initialize the structures. Each structure will have as many columns
    % as there are genotypes. Cell arrays within each coulm will have size
    % of MxN where M is the number of flies and N is number of total
    % conditions (dark + 8 different light conditions as of 220131).
    % Sequence of light conditions will be stored in the seq variable below
    lightCondNum = 2;
    VfwThrS = 14; % mm/s | the maximum forward velocity that will be
    % considered as slow walking
    VfwThrF = 18; % mm/s | the minimum forward velocity that will be
    % considered as fast walking
    AD = struct();
    AD.VfwThr = [VfwThrS,VfwThrF];
    % structure fields for all forward bouts
    AD.Wmu = cell(1,length(genotypes)); % weighted mean
    AD.Wstd = cell(1,length(genotypes)); % std of the weighted mean
    AD.AbsWmu = cell(1,length(genotypes)); % weighted mean of abs angdev
    AD.BWmu = cell(1,length(genotypes)); %weighted mean of binarized angdev
    AD.BlSum = cell(1,length(genotypes)); % bout length summed per fly
    AD.Bl = cell(1,length(genotypes)); % bout length per bout per fly
    AD.Adev = cell(1,length(genotypes)); % bout length per bout per fly
    AD.Blmu = cell(1,length(genotypes)); % mean bout length per fly
    AD.Vfmu = cell(1,length(genotypes)); % mean Vfw per bout per fly
    AD.VFperGen = cell(1,length(genotypes)); % Vfw values per genotype
    % structure fields for "fast" forward bouts
    AD.WmuF = cell(1,length(genotypes)); % weighted mean
    AD.WstdF = cell(1,length(genotypes)); % std of the weighted mean
    AD.AbsWmuF = cell(1,length(genotypes)); % weighted mean of abs angdev
    AD.BlSumF = cell(1,length(genotypes)); % bout length summed per fly
    AD.BlF = cell(1,length(genotypes)); % bout length per bout per fly
    AD.BlmuF = cell(1,length(genotypes)); % mean bout length per fly
    AD.VfmuF = cell(1,length(genotypes)); % mean Vfw per bout per fly
    % structure fields for "slow" forward bouts
    AD.WmuS = cell(1,length(genotypes)); % weighted mean
    AD.WstdS = cell(1,length(genotypes)); % std of the weighted mean
    AD.AbsWmuS = cell(1,length(genotypes)); % weighted mean of abs angdev
    AD.BlSumS = cell(1,length(genotypes)); % bout length summed per fly
    AD.BlS = cell(1,length(genotypes)); % bout length per bout per fly
    AD.BlmuS = cell(1,length(genotypes)); % mean bout length per fly
    AD.VfmuS = cell(1,length(genotypes)); % mean Vfw per bout per fly
    % All angular deviation values concatenated per genotype
    AD.AdevAll = cell(1,length(genotypes));

    for c = 1:length(cond)
        for g = 1:length(genotypes) % for each genotype
        % directory where the files are located
        path = [baseDir,behDir,cond{c},'\',genotypes{g},'\'];
        disp([behDir,cond{c},'\',genotypes{g},'\'])
        % load the low resolution data relative to the forward segments
        % together with the stimuli sequence
        [VForwSeg,seq] = GetForwSegVData(path, params);
        % NOTE: Fw speed cutoff is stored in params.vft and it can be
        % changed based on slow/fast forward speed analysis preference

        % binning windows for the angular displacements (deg/mm)
        DevCents = [-1000 -3.75 -2.5 -1.25 0 1.25 2.5 3.75 1000];

        % initialize variables
        VR = cell(length(DevCents)-1,1); % rotational velocity
        VF = cell(length(DevCents)-1,1); % forward velocity
        VS = cell(length(DevCents)-1,1); % side velocity
        NB = cell(length(DevCents)-1,1); % number of bouts
        angDevAll = [];
        % Length of each forward walking bout per fly
        Blength = cell(size(VForwSeg,2),size(VForwSeg,1));
        % Angular deviation per bout per fly
        Adev = cell(size(VForwSeg,2),size(VForwSeg,1)); % signed
        AdevAbs = cell(size(VForwSeg,2),size(VForwSeg,1)); % absolute
        AdevB = cell(size(VForwSeg,2),size(VForwSeg,1)); % binarized
        % Mean forward velocity per bout per fly
        Vfw = cell(size(VForwSeg,2),size(VForwSeg,1));

        for jk = 1:length(DevCents)-1 % for each angular deviation bin
            % define boundaries of each bin
            minAngDev = DevCents(jk);
            maxAngDev = DevCents(jk+1);

            for s = 1:size(VForwSeg,1) % for each stimulus
                for n = 1:size(VForwSeg,2) % for each fly
                    for i = 1:length(VForwSeg{s,n}) % for each fw bout
                        if ~isempty(VForwSeg{s,n}{i})% if there is data
                            % low res speed data for this fw segment
                            vr = VForwSeg{s,n}{i}.VrLR; %Vrotational
                            vf = VForwSeg{s,n}{i}.VfLR; %Vforward
                            vs = VForwSeg{s,n}{i}.VsLR; %Vsideways
                            vt = sqrt(vf.*vf + vs.*vs); %Vtranslational
                            angDev = mean(vr./vt); % angular deviation
                            angDevB = angDev/abs(angDev); % binarized

                            %if the angular deviation belong to the bin
                            if angDev > minAngDev && angDev < maxAngDev
                                % store the data appropriately
                                VR{jk} = vertcat(VR{jk}, vr);
                                VF{jk} = vertcat(VF{jk}, vf);
                                VS{jk} = vertcat(VS{jk}, vs);
                                NB{jk} = vertcat(NB{jk}, length(vs));
                                angDevAll = vertcat(angDevAll,angDev);
                                Adev{n,s} = vertcat(Adev{n,s},angDev);
                                AdevB{n,s} = vertcat(AdevB{n,s},angDevB);
                                AdevAbs{n,s} = vertcat(AdevAbs{n,s},...
                                    abs(angDev));
                                Vfw{n,s} = vertcat(Vfw{n,s},mean(vf));
                                Blength{n,s} = vertcat(Blength{n,s},...
                                    length(vs));
                            end
                        end
                    end
                end
            end
        end

        % initialize variables
        % Angular deviation per fly weighted by ALL fw bout length
        AngDevWmu = nan(size(Adev));
        AngDevWvar = nan(size(Adev));
        AngDevAbsWmu = nan(size(Adev));
        AngDevBWmu = nan(size(AdevB));
        BlSum = nan(size(Adev));
        BlW = nan(size(Adev));
        AD.AdevAll{g} = vertcat(AD.AdevAll{g},angDevAll);
        % Angular deviation per fly weighted by FAST fw bout length
        AngDevWmuF = nan(size(Adev));
        AngDevWvarF = nan(size(Adev));
        AngDevAbsWmuF = nan(size(Adev));
        BlSumF = nan(size(Adev));
        AdevF = cell(size(Adev));
        AdevAbsF = cell(size(AdevAbs));
        BlengthF = cell(size(Blength));
        % Angular deviation per fly weighted by SLOW fw bout length
        AngDevWmuS = nan(size(Adev));
        AngDevWvarS = nan(size(Adev));
        AngDevAbsWmuS = nan(size(Adev));
        BlSumS = nan(size(Adev));
        AdevS = cell(size(Adev));
        BlengthS = cell(size(Blength));
        AdevAbsS = cell(size(AdevAbs));

        % logical index arrays for slow and fast bouts
        iFast =cellfun(@(x) x >= VfwThrF, Vfw, 'UniformOutput', false);
        iSlow =cellfun(@(x) x < VfwThrS, Vfw, 'UniformOutput', false);

        for n = 1:size(Adev,1) % for each fly
            for i = 1:size(Adev,2) % for each stimulus

                % Angular deviations of all bouts calculated as a
                % weighted mean based on bout lengths
                AngDevWmu(n,i)=...
                    nansum((Adev{n,i}.*Blength{n,i})/...
                    nansum(Blength{n,i}));
                AngDevAbsWmu(n,i)=...
                    nansum((AdevAbs{n,i}.*Blength{n,i})/...
                    nansum(Blength{n,i}));
                AngDevBWmu(n,i)=...
                    nansum((AdevB{n,i}.*Blength{n,i})/...
                    nansum(Blength{n,i}));
                AngDevWvar(n,i)=...
                    nansum((((Adev{n,i}-AngDevWmu(n,i)).^2)...
                    .*Blength{n,i})/nansum(Blength{n,i}));
                BlSum(n,i) = nansum(Blength{n,i});

                % Angular deviations of fast bouts calculated as a
                % weighted mean based on bout lengths
                AdevF{n,i} = Adev{n,i}(iFast{n,i}); % fast bouts only
                AdevAbsF{n,i} = AdevAbs{n,i}(iFast{n,i});
                BlengthF{n,i} = Blength{n,i}(iFast{n,i});
                AngDevWmuF(n,i)=...
                    nansum((AdevF{n,i}.*BlengthF{n,i})/...
                    nansum(BlengthF{n,i}));
                AngDevAbsWmuF(n,i)=...
                    nansum((AdevAbsF{n,i}.*BlengthF{n,i})/...
                    nansum(BlengthF{n,i}));
                AngDevWvarF(n,i)=...
                    nansum((((AdevF{n,i}-AngDevWmuF(n,i)).^2)...
                    .*BlengthF{n,i})/nansum(BlengthF{n,i}));
                BlSumF(n,i) = nansum(BlengthF{n,i});

                % Angular deviations of slow bouts calculated as a
                % weighted mean based on bout lengths
                AdevS{n,i} = Adev{n,i}(iSlow{n,i}); % slow bouts only
                AdevAbsS{n,i} = AdevAbs{n,i}(iSlow{n,i});
                BlengthS{n,i} = Blength{n,i}(iSlow{n,i});
                AngDevWmuS(n,i)=...
                    nansum((AdevS{n,i}.*BlengthS{n,i})/...
                    nansum(BlengthS{n,i}));
                AngDevAbsWmuS(n,i)=...
                    nansum((AdevAbsS{n,i}.*BlengthS{n,i})/...
                    nansum(BlengthS{n,i}));
                AngDevWvarS(n,i)=...
                    nansum((((AdevS{n,i}-AngDevWmuS(n,i)).^2)...
                    .*BlengthS{n,i})/nansum(BlengthS{n,i}));
                BlSumS(n,i) = nansum(BlengthS{n,i});
            end
        end

        % remove zeros from the data
        AngDevWmu(BlSum==0) = nan;
        AngDevWvar(BlSum==0) = nan;
        AngDevAbsWmu(BlSum==0) = nan;
        AngDevBWmu(BlSum==0) = nan;
        BlSum(BlSum==0) = nan;
        AngDevWmuF(BlSumF==0) = nan;
        AngDevWvarF(BlSumF==0) = nan;
        AngDevAbsWmuF(BlSumF==0) = nan;
        BlSumF(BlSumF==0) = nan;
        AngDevWmuS(BlSumS==0) = nan;
        AngDevWvarS(BlSumS==0) = nan;
        AngDevAbsWmuS(BlSumS==0) = nan;
        BlSumS(BlSumS==0) = nan;
        % standard deviation of the weighted means
        AngDevWstd = sqrt(AngDevWvar);
        AngDevWstdF = sqrt(AngDevWvarF);
        AngDevWstdS = sqrt(AngDevWvarS);

        % horzcat concatenates in the opposite of what we want when the
        % initial array is empty, so we swap the positions of variables
        % during the first concatenation (which happens in 'Dark'
        % condition). If you get errors here, make sure that the first
        % element of conds variable is 'Dark'.
        switch cond{c}
            case 'Dark'
                AD.Wmu{g} = horzcat(AngDevWmu,AD.Wmu{g});
                AD.Wstd{g} = horzcat(AngDevWstd,AD.Wstd{g});
                AD.AbsWmu{g} = horzcat(AngDevAbsWmu,AD.AbsWmu{g});
                AD.BWmu{g} = horzcat(AngDevBWmu,AD.BWmu{g});
                AD.BlSum{g} = horzcat(BlSum,AD.BlSum{g});
                AD.Bl{g} = horzcat(Blength,AD.Bl{g});
                AD.Adev{g} = horzcat(Adev,AD.Adev{g});
                AD.WmuF{g} = horzcat(AngDevWmuF,AD.WmuF{g});
                AD.WstdF{g} = horzcat(AngDevWstdF,AD.WstdF{g});
                AD.AbsWmuF{g} = horzcat(AngDevAbsWmuF,AD.AbsWmuF{g});
                AD.BlSumF{g} = horzcat(BlSumF,AD.BlSumF{g});
                AD.BlF{g} = horzcat(BlengthF,AD.BlF{g});
                AD.WmuS{g} = horzcat(AngDevWmuS,AD.WmuS{g});
                AD.WstdS{g} = horzcat(AngDevWstdS,AD.WstdS{g});
                AD.AbsWmuS{g} = horzcat(AngDevAbsWmuS,AD.AbsWmuS{g});
                AD.BlSumS{g} = horzcat(BlSumS,AD.BlSumS{g});
                AD.BlS{g} = horzcat(BlengthS,AD.BlS{g});

            case {'Light','Light Bilateral'}
                AD.Wmu{g} = horzcat(AD.Wmu{g},AngDevWmu);
                AD.Wstd{g} = horzcat(AD.Wstd{g},AngDevWstd);
                AD.AbsWmu{g} = horzcat(AD.AbsWmu{g},AngDevAbsWmu);
                AD.BWmu{g} = horzcat(AD.BWmu{g},AngDevBWmu);
                AD.BlSum{g} = horzcat(AD.BlSum{g},BlSum);
                AD.Bl{g} = horzcat(AD.Bl{g},Blength);
                AD.Adev{g} = horzcat(AD.Adev{g},Adev);
                AD.WmuF{g} = horzcat(AD.WmuF{g},AngDevWmuF);
                AD.WstdF{g} = horzcat(AD.WstdF{g},AngDevWstdF);
                AD.AbsWmuF{g} = horzcat(AD.AbsWmuF{g},AngDevAbsWmuF);
                AD.BlSumF{g} = horzcat(AD.BlSumF{g},BlSumF);
                AD.BlF{g} = horzcat(AD.BlF{g},BlengthF);
                AD.WmuS{g} = horzcat(AD.WmuS{g},AngDevWmuS);
                AD.WstdS{g} = horzcat(AD.WstdS{g},AngDevWstdS);
                AD.AbsWmuS{g} = horzcat(AD.AbsWmuS{g},AngDevAbsWmuS);
                AD.BlSumS{g} = horzcat(AD.BlSumS{g},BlSumS);
                AD.BlS{g} = horzcat(AD.BlS{g},BlengthS);
                AD.seq = horzcat({'Dark'},seq'); % add stimulus info
        end
        
        %% Extended Data Figure 10a - Plot example trajectories
        xT = cell(1,1);
        yT = cell(1,1);
        indd = cell(1,1);
        a = 1;
        cmap1 = jet(length(VR));
        xTall = [];
        yTall = [];
        % iterate through the angular deviation bins
        for i = 2 : length(NB)-1
            nbs = NB{i};
            nbScs = [1; cumsum(nbs)];
            % iterate through the forward segments
            for j = 1 : length(nbs)
                vr = VR{i}(nbScs(j):nbScs(j+1))/60;
                vf = VF{i}(nbScs(j):nbScs(j+1))/60;
                vs = VS{i}(nbScs(j):nbScs(j+1))/60;
                x = [];
                y = [];
                th = [];
                if length(vr) < 180
                    % get centered X and Y
                    for ij = 1 : length(vr)
                        if isempty(x)
                            th = vertcat(th,pi*vr(ij)/180);
                            x = vertcat(x,vf(ij)*sin(0)+vs(ij)*cos(0));
                            y = vertcat(y,vf(ij)*cos(0)-vs(ij)*sin(0));
                        else
                            th = vertcat(th,th(end)+pi*vr(ij)/180);
                            x = vertcat(x,x(end)+...
                                vf(ij)*sin(th(ij-1))...
                                +vs(ij)*cos(th(ij-1)));
                            y = vertcat(y,y(end)+...
                                vf(ij)*cos(th(ij-1))...
                                -vs(ij)*sin(th(ij-1)));
                        end
                    end
                end
                xT{a} = x;
                yT{a} = y;
                indd{a} = i;
                if a == 1
                    xTall = x';
                    yTall = y';
                else
                    if length(x) > size(xTall,2)
                        xTall(:,(size(xTall,2)+1):length(x)) = NaN;
                        yTall(:,(size(yTall,2)+1):length(y)) = NaN;
                        xTall(a,:) = x';
                        yTall(a,:) = y';
                    else
                        xTall(a,:) = NaN;
                        yTall(a,:) = NaN;
                        xTall(a,1:length(x)) = x';
                        yTall(a,1:length(y)) = y';
                    end
                end
                a = a + 1;
            end
        end
        yTallBin = 0:round(max(nanmax(yTall)));
        xTallBin = nan(1,length(yTallBin));
        yTallR = round(yTall);
        for i = 1:(round(max(nanmax(yTall)))+1)
            xTallBin(i)=nanmean(xTall(yTallR==i-1));
        end
        % plot a color coded sample of the forward segments
        inds = randperm(length(xT));
        frac = 100; % percentage of data to be plotted
        ythr = 10; % threshold of forward translation (mm in y axis)
        inds = inds(1:round(length(xT)*frac/100));

        fs = 18; % font size
        f=figure('WindowState','maximized');
        hold on
        for i = 1 : length(inds)
            if max(yT{inds(i)})>ythr
                plot(yT{inds(i)}, xT{inds(i)},...
                    'color', cmap1(indd{inds(i)},:),'LineWidth',1.5)
            end
        end
        plot(yTallBin,xTallBin,'color','k','LineWidth',2)

        axis([-2 30 -25 25])
        limsX = xlim;
        line(limsX,[0 0],'Color','k','LineWidth',0.5,'LineStyle','--')
        %             title(['bIPS-',gTypes{g}]) % add title, tick and axis labels
        %             title(['DNp15-',genotypes{g}]) % add title, tick and axis labels
        title(['DNp15-RicinA-',gTypes{g}]) % add title, tick and axis labels
        ylabel('X position (mm)'); xlabel('Y position (mm)')
        set(gca,'FontSize',fs)
        caxis([DevCents(2) DevCents(end-1)])
        colormap(cmap1)
        ch = colorbar('Ticks',DevCents(2:end-1));
        ch.Label.String = 'Ang. Dev. (deg/mm)';
        ch.Label.Rotation = 270;
        ts = timestamp(); % get current date
        if SaveFigures
            set(gcf,'renderer','Painters')
            print(f,strcat(saveDir,ts,'-Paths-',...
                genotypes{g},'-',cond{c},'-',num2str(frac),...
                'perc-yThr',num2str(ythr)),'-dtiff') %#ok<*UNRCH>
            print(f,strcat(saveDir,ts,'-Paths-',...
                genotypes{g},'-',cond{c},'-',num2str(frac),...
                'perc-yThr',num2str(ythr)),'-depsc')
            close(f)
        end
        end
    end
    
    %     AD.AdevAll = AdevAll;
    for g = 1:length(genotypes) % for each genotype
        nFly = size(AD.Bl{1,g}(:,1),1);
        Blmean = nan(nFly,size(AD.seq,2));
        BlmeanF = nan(nFly,size(AD.seq,2));
        BlmeanS = nan(nFly,size(AD.seq,2));
        for i = 1:size(AD.seq,2) % for each condition
            for n = 1:nFly % for each fly
                Blmean(n,i) = nanmean(AD.Bl{1,g}{n,i});
                BlmeanF(n,i) = nanmean(AD.BlF{1,g}{n,i});
                BlmeanS(n,i) = nanmean(AD.BlS{1,g}{n,i});
            end
        end
        AD.Blmu{1,g} = Blmean;
        AD.BlmuF{1,g} = BlmeanF;
        AD.BlmuS{1,g} = BlmeanS;
    end
    if ~exist(saveDir, 'dir')
        mkdir(saveDir);
    end
    save([saveDir,'AngDevVars.mat'],'AD')
else
    load([saveDir,'AngDevVars.mat'])
end


% Plot weigthed mean angular deviation variables per condition
% --------------------------- User Input ----------------------------------
% Here you can specify which angular deviation related data you'd like to
% plot for visualization. All data weighted by total bout length per fly.
% 'Mean' will plot the mean angular deviation per fly
% 'Abs' will plot the mean absolute angular deviation per fly
% 'Binary' will plot the binarized angular deviation per fly
% 'Std' will plot the standard deviation of mean angular deviation per fly
% 'Bout' will plot the mean forward bout lengths per fly
% Accepted arguments: ['Mean','Abs','Binary','Std','Bout']
AngDevToPlot = 'Mean';

% Here you can filter the data plotted based on the mean forward velocity
% of the fly during a forward bout. As of 221130, definitions of fast and
% slow bouts are static and pre-determined during data loading phase.
% 'Slow' and 'Fast' will plot only the slow and fast bouts, 'All' will plot
% all the data. Note that Slow + Fast does not have to equal All bouts.
% Accepted arguments: ['Slow','Fast','All']
SpeedToPlot = 'All';

sig = 1; % 1 if you want to overlay significance stars (MWW non-param test)
fs = 14; % font size for plot labels
alpha = 0.5; % transparency of the dots
ms = 32; % dot size of the mean data point
scaleFactor = 0.1; % scale the dot size
scatterFactor = 0.26; % add scatter (fraction between 0 and 1) to the dots
limsX = [0 size(AD.Wmu,2)+1]; % axis limits

%----------------------------end of user input-----------------------------

switch size(AD.seq,2) % change plotting order so 10deg and RG comes last
    case 1 % if there is only dark
        conds = 1; xTicks = 1;
    case 2
        % 1 + 10 NG  (bIPS Translational)
        conds = [2,1]; xTicks = 1:size(AD.Wmu,2);
    case 3
        % Dark + 1 + 10 NG (DNp15 Ricin A)
        conds = [1,3,2]; xTicks = 1:size(AD.Wmu,2);
    case 7        
        % Dark + 1 + 5 + 10 NG|RG  (DNp15 hsFLP)
        conds = [1,4,6,2,5,7,3]; xTicks = 1:size(AD.Wmu,2);
    case 8
        conds = [3,5,7,1,4,6,8,2]; xTicks = 1:8;
    case 9 % assumes dark is present
        % Dark + 1 + 5 + 10 NG|RG  (bIPS hsFLP)
        % conds = [1,4,6,8,2,5,7,9,3]; xTicks = 1:size(AD.Wmu,2);
        conds = [1,4,2]; xTicks = 1:3;     
end

% calculate maximum y values across data points to have the significance
% star plots well separated over the data points
my = zeros(1,length(conds));

f=figure('WindowState','maximized');

i = find(contains(AD.seq,'10NG'));
sf = scaleFactor;
hold on;

if contains(behDir,'DNp15')
    col = cReg(1,:);
elseif contains(behDir,'bIPS')
    col = cReg(2,:);
else
    col = [0,0,0];
end

for g = 1:length(genotypes)

    switch AngDevToPlot
        % vp: variable to plot | yl: y-axis label | sn: save name |
        case 'Mean'
            switch SpeedToPlot
                case 'Slow'
                    vp = AD.WmuS;
                case 'Fast'
                    vp = AD.WmuF;
                case 'All'
                    vp = AD.Wmu;
                otherwise
                    error('Unknown/unspecified SpeedToPlot variable!')
            end
            sn = 'AngDevPerCond';
            yl = 'Ang. Dev. (deg/mm)';
            limsY = [-7 7];
        case 'Abs'
            switch SpeedToPlot
                case 'Slow'
                    vp = AD.AbsWmuS;
                case 'Fast'
                    vp = AD.AbsWmuF;
                case 'All'
                    vp = AD.AbsWmu;
                otherwise
                    error('Unknown/unspecified SpeedToPlot variable!')
            end
            sn = 'AbsAngDevPerCond';
            yl = 'Abs. Ang. Dev. (deg/mm)';
            limsY = [0 5];
        case 'Std'
            switch SpeedToPlot
                case 'Slow'
                    vp = AD.WstdS;
                case 'Fast'
                    vp = AD.WstdF;
                case 'All'
                    vp = AD.Wstd;
                otherwise
                    error('Unknown/unspecified SpeedToPlot variable!')
            end
            sn = 'StdAngDevPerCond';
            yl = 'std of Ang. Dev. (a.u.)';
            limsY = [0 10];
        case 'Bout'
            switch SpeedToPlot
                case 'Slow'
                    vp = AD.BlmuS;
                case 'Fast'
                    vp = AD.BlmuF;
                case 'All'
                    vp = AD.Blmu;
                otherwise
                    error('Unknown/unspecified SpeedToPlot variable!')
            end
            sn = 'FwBoutLengthPerCond';
            yl = 'Fw bout length (frames)';
            limsY = [20 80];
        case 'Binary'
            switch SpeedToPlot
                case 'Slow'
                    error('Binarized plots for Slow not implemented!')
                case 'Fast'
                    error('Binarized plots for Fast not implemented!')
                case 'All'
                    vp = AD.BWmu;
                otherwise
                    error('Unknown/unspecified SpeedToPlot variable!')
            end
            sn = 'AngBiasPerCond';
            yl = 'Angular Bias';
            limsY = [-1.1 1.1];
        otherwise
            error(['Don''t know what to plot, ',...
                'please check the AngDevToPlot variable!'])
    end
    switch SpeedToPlot
        % bl: bout length data used to weigh the per fly data
        case 'Slow'
            bl = AD.BlSumS;
        case 'Fast'
            bl = AD.BlSumF;
        case 'All'
            bl = AD.BlSum;
        otherwise
            error('Unknown/unspecified SpeedToPlot variable!')
    end

    % calculate weighted grand mean and its variance & std deviation
    Wmean = nansum((vp{1,gi(g)}(:,i).*...
        bl{1,gi(g)}(:,i))/nansum(bl{1,gi(g)}(:,i)));
    Wvar = nansum((((vp{1,gi(g)}(:,i) - Wmean).^2).*...
        bl{1,gi(g)} (:,i))/nansum(bl{1,gi(g)}(:,i)));
    Wstd = sqrt(Wvar);
    % define x-axis location of each dot
    sc = (rand(size(AD.Blmu{1,gi(g)},1),1)-0.5)*scatterFactor;
    x = (ones(size(AD.Blmu{1,gi(g)},1),1)*g) + sc;
    % plot the data, add mean and std error bars
    scatter(x,vp{1,gi(g)}(:,i),...
        (bl{1,gi(g)}(:,i)*sf),...
        'MarkerFaceColor',col,'MarkerEdgeColor','none',...
        'MarkerFaceAlpha',alpha);
    plot(g,Wmean,'k.','MarkerSize',ms)
    nf = size(AD.Blmu{1,gi(g)},1); % number of flies in this genotype
    eh = errorbar(g,Wmean,Wstd/sqrt(nf));
    set(eh,'LineWidth',2,'Color','k')
    maxbl = nanmax(vp{1,gi(g)}(:,i));
    my(i) = max([my(i),maxbl,limsY(2)*0.4]);
    my(i) = min(my(i),limsY(2)*0.4);
end
if sig
    % Perform Mann-Whitney-Wilcoxon non parametric test and plot
    % significance indicators if p<0.05
    nck = nchoosek(1:g,2);
    %         % bIPS TrpA1 Kir combined plots only: remove unnecessary comparisons
    %         nck = nck([1,10:end],:);
    nph = size(nck,1); % number of comparisons made (used for adjusting the
    % minimum P value required for significance - a.k.a. Bonferroni correction)
    for ii = 1 : size(nck,1)
        pair1 = vp{1,gi(nck(ii,1))}(:,i);
        pair2 = vp{1,gi(nck(ii,2))}(:,i);
        pair1(isnan(pair1)) = []; pair2(isnan(pair2)) = [];
        stats{ii} = mwwtest(pair1,pair2,0);
        % uncomment below if you want to display all p-values per pair
        % disp([AD.seq{i},'-',gTypes{nck(ii,1)},'-',gTypes{nck(ii,2)},'-',num2str(stats{ii}.p(2))])
        if stats{ii}.p < 0.05/nph
            % uncomment below if you want to display significant p-values per pair
            %                 disp([AD.seq{i},'-',gTypes{nck(ii,1)},'-',gTypes{nck(ii,2)},'-',num2str(stats{ii}.p(2))])
            sigline([(nck(ii,1)),(nck(ii,2))],num2str(stats{ii}.p(2)),my(i)+my(i)*0.1*ii);
        end
    end
end
% formatting the plot
% set axis limits
xlim(limsX); ylim(limsY);
% add tick and axis labels
ylabel(yl)
set(gca,'XTick',xTicks,'XTickLabel',gOrder,'FontSize',fs,...
    'XTickLabelRotation',45)
if contains(AngDevToPlot,'Mean') || contains(AngDevToPlot,'Binary')
    % add dashed line at 0
    line(limsX,[0 0],'Color','k','LineWidth',0.5,'LineStyle','--')
end
% text overlay of speed
xL = 2; % text location on the x-axis
yL = my(i)-my(i)*0.11*ii; % text location on the y-axis
switch SpeedToPlot
    case 'Slow'
        T = text(xL,yL,['Slow (<',num2str(AD.VfwThr(1)),'mm/s)']);
    case 'Fast'
        T = text(xL,yL,['Fast (>',num2str(AD.VfwThr(2)),'mm/s)']);
end

%% Extended Data Figure 10c
% DNp15 hsFlp RicinA
gTypes = {'Both','Left','None','Right'};
behDir = 'DNp15Ricin\';
saveDir = [baseDir,'Analysis\DNp15Ricin\']; % save directory
path = [baseDir,behDir,'Light'];
[~,gi] = ismember(gOrder,gTypes); % plotting order indices

if exist([saveDir,'SaccadeVars11.mat'],'file')
    load([saveDir,'SaccadeVars11.mat'])
else
    cd(path);
    datedir = dir;
    filenames = {datedir(3:end).name};
    for k = 1 : length(filenames)
        pathF{k} = [path,'\',filenames{k},'\'];

        if ~isempty(pathF{k})
            [SPKVR{k}, ~, SPKVRDist{k},SPKVF{k}, pcav{k}, pTypes{k}, dt{k}] = GetSpikePCA(pathF{k}, params);
        end

        % visual stimuli order. Useful to plot in numerical order
        % (e.g. 1NG 5NG 10NG instead of 10NG 1NG 5NG)
        if strcmp('Dark',genotypes{k})
            seq{k} = 1;
        else
            seq{k} = [2,1];
        end
    end
end

% Plot saccade direction bias
fs = 14; % font size for plot labels
alpha = 0.7; % transparency of the dots
ms=32; % dot size of the mean
limsX = [0.3 length(genotypes)+0.7]; % axis limits
limsY = [-0.7 0.7];
% limsY = [-1 1];
scatterFactor = 0.16; % add scatter (fraction between 0 and 1) to the dots

switch cond{1}
    case 'Dark'
        xTicks = 1:length(genotypes);
        colors = [22,0,22]/256; conds = 1;
    case {'Light','Light Bilateral'}
        conds = seq{1};
        xTicks = 1:length(genotypes);
end

if contains(behDir,'DNp15')
    col = cReg(1,:);
elseif contains(behDir,'bIPS')
    col = cReg(2,:);
else
    col = [0,0,0];
end

my = zeros(1,length(conds));
SPKbinary = cell(size(SPKVR));
SPKnum = cell(size(SPKVR));
f=figure('WindowState','maximized');

c = find(contains(pTypes{1},'1NG'));
hold on;
for g = 1:length(genotypes) % for each genotype
    nFly = size(SPKVR{1,gi(g)},1);
    SPKbinary{1,g} = nan(nFly,length(cond(1)));
    for fly = 1:nFly % for each fly
        spikes = SPKVR{1,gi(g)}{fly,seq{1}(c)};
        sbinary = nan(1,size(spikes,2));
        sbinary(max(spikes)<max(abs(spikes))) = -1;
        sbinary(max(spikes)==max(abs(spikes))) = 1;
        sbinary(spikes>0) = 1;
        sbinary(spikes<0) = -1;
        SPKbinary{1,g}(fly,seq{1}(c)) = sum(sbinary)/length(sbinary);
        SPKnum{1,g}(fly,seq{1}(c)) = length(sbinary);
    end
    vp = SPKbinary{1,g}(:,seq{1}(c)); % value to plot
    vpn = SPKnum{1,g}(:,seq{1}(c)); % num of events (for weighted vp)
    % remove flies without any saccades
    nd = sum(vpn==0); % number of discarded flies
    vp(vpn==0) = []; vpn(vpn==0) = [];

    Wmean = nansum((vp.*vpn)/sum(vpn));
    Wvar = nansum((((vp - Wmean).^2).*vpn)/sum(vpn));
    Wstd = sqrt(Wvar);

    sc = (rand(size(vp,1),1)-0.5)*scatterFactor;
    x = (ones(size(vp,1),1)*g) + sc;
    % plot the data, add mean and std error bars
    scatter(x,vp,vpn,...
        'MarkerFaceColor',col,'MarkerEdgeColor','none',...
        'MarkerFaceAlpha',alpha);

    plot(g,Wmean,'k.','MarkerSize',ms)
    eh=errorbar(g,Wmean,(1.960*Wstd/sqrt(length(vp))));
    set(eh,'LineWidth',2,'Color','k')

    maxbl = nanmax(SPKbinary{1,g}(:,seq{1}(c)));
    my(seq{1}(c)) = max(my(seq{1}(c)),maxbl);
end
% set axis limits
xlim(limsX);
ylim(limsY);
ylabel('Saccade bias')
set(gca,'XTick',xTicks,'XTickLabel',gOrder,'FontSize',fs,...
    'XTickLabelRotation',45)
% add dashed line at 0
line(limsX,[0 0],'Color','k','LineWidth',0.5,'LineStyle','--')

nck = nchoosek(1:g,2);
nph = size(nck,1); % number of comparisons made (used for adjusting the
% minimum P value required for significance - a.k.a. Bonferroni correction)
for ii = 1 : size(nck,1)
    pair1 = SPKbinary{1,nck(ii,1)}(:,seq{1}(c));
    pair2 = SPKbinary{1,nck(ii,2)}(:,seq{1}(c));
    pair1(isnan(pair1)) = []; pair2(isnan(pair2)) = [];
    stats{ii} = mwwtest(pair1,pair2,0);
    if stats{ii}.p(2) < 0.05/nph
        sigline([nck(ii,1),nck(ii,2)],num2str(stats{ii}.p(2)),my(seq{1}(c))+my(seq{1}(c))*0.05*ii);
    end
end

%% Extended Data Figure 10d and 10f
% DNp15 hsFlp Kir
cond = {'Dark','Light'}; % experimental conditions
genotypes = {'No cells','Single cell Left','Single Cell Right','Double cell'};
gTypes = {'None','Left','Right','Both'};
behDir = 'DNp15Kir\';
saveDir = [baseDir,'Analysis\DNp15Kir\']; % save directory
path = [baseDir,behDir,'Light'];
controls = [0 0 1 0];

[~,gi] = ismember(gOrder,gTypes); % plotting order indices


if exist([saveDir,'SaccadeVars11.mat'],'file')
    load([saveDir,'SaccadeVars11.mat'])
    [~,gi] = ismember(gOrder,gTypes); % plotting order indices
else
    cd(path);
    datedir = dir;
    filenames = {datedir(3:end).name};
    for k = 1 : length(filenames)
        pathF{k} = [path,'\',filenames{k},'\'];
        if ~isempty(pathF{k})
            [SPKVR{k}, ~, SPKVRDist{k},SPKVF{k}, pcav{k}, pTypes{k}, dt{k}] = GetSpikePCA(pathF{k}, params);
        end
        % visual stimuli order. Useful to plot in numerical order
        % (e.g. 1NG 5NG 10NG instead of 10NG 1NG 5NG)
        if strcmp('Dark',genotypes{k})
            seq{k} = 1;
        else
            seq{k} = [2,1];
        end
    end
end

% Extended Data Figure 10d - Plot saccade direction bias
fs = 14; % font size for plot labels
alpha = 0.7; % transparency of the dots
ms=32; % dot size of the mean
limsX = [0.3 length(genotypes)+0.7]; % axis limits
limsY = [-0.7 0.7];
% limsY = [-1 1];
scatterFactor = 0.16; % add scatter (fraction between 0 and 1) to the dots

switch cond{1}
    case 'Dark'
        xTicks = 1:length(genotypes);
        colors = [22,0,22]/256; conds = 1;
    case {'Light','Light Bilateral'}
        conds = seq{1};
        xTicks = 1:length(genotypes);
end

if contains(behDir,'DNp15')
    col = cReg(1,:);
elseif contains(behDir,'bIPS')
    col = cReg(2,:);
else
    col = [0,0,0];
end

my = zeros(1,length(conds));
SPKbinary = cell(size(SPKVR));
SPKnum = cell(size(SPKVR));
f=figure('WindowState','maximized');

c = find(contains(pTypes{1},'10NG'));
hold on;
for g = 1:length(genotypes) % for each genotype
    nFly = size(SPKVR{1,gi(g)},1);
    SPKbinary{1,g} = nan(nFly,length(cond(1)));
    for fly = 1:nFly % for each fly
        spikes = SPKVR{1,gi(g)}{fly,c};
        sbinary = nan(1,size(spikes,2));
        sbinary(max(spikes)<max(abs(spikes))) = -1;
        sbinary(max(spikes)==max(abs(spikes))) = 1;
        sbinary(spikes>0) = 1;
        sbinary(spikes<0) = -1;
        SPKbinary{1,g}(fly,c) = sum(sbinary)/length(sbinary);
        SPKnum{1,g}(fly,c) = length(sbinary);
    end
    vp = SPKbinary{1,g}(:,c); % value to plot
    vpn = SPKnum{1,g}(:,c); % num of events (for weighted vp)
    % remove flies without any saccades
    nd = sum(vpn==0); % number of discarded flies
    vp(vpn==0) = []; vpn(vpn==0) = [];

    Wmean = nansum((vp.*vpn)/sum(vpn));
    Wvar = nansum((((vp - Wmean).^2).*vpn)/sum(vpn));
    Wstd = sqrt(Wvar);

    sc = (rand(size(vp,1),1)-0.5)*scatterFactor;
    x = (ones(size(vp,1),1)*g) + sc;
    % plot the data, add mean and std error bars
    scatter(x,vp,vpn,...
        'MarkerFaceColor',col,'MarkerEdgeColor','none',...
        'MarkerFaceAlpha',alpha);

    plot(g,Wmean,'k.','MarkerSize',ms)
    eh=errorbar(g,Wmean,(1.960*Wstd/sqrt(length(vp))));
    set(eh,'LineWidth',2,'Color','k')

    maxbl = nanmax(SPKbinary{1,g}(:,c));
    my(c) = max(my(c),maxbl);
end
% set axis limits
xlim(limsX);
ylim(limsY);
ylabel('Saccade bias')
set(gca,'XTick',xTicks,'XTickLabel',gOrder,'FontSize',fs,...
    'XTickLabelRotation',45)
% add dashed line at 0
line(limsX,[0 0],'Color','k','LineWidth',0.5,'LineStyle','--')

nck = nchoosek(1:g,2);
nph = size(nck,1); % number of comparisons made (used for adjusting the
% minimum P value required for significance - a.k.a. Bonferroni correction)
for ii = 1 : size(nck,1)
    pair1 = SPKbinary{1,nck(ii,1)}(:,c);
    pair2 = SPKbinary{1,nck(ii,2)}(:,c);
    pair1(isnan(pair1)) = []; pair2(isnan(pair2)) = [];
    stats{ii} = mwwtest(pair1,pair2,0);
    if stats{ii}.p(2) < 0.05/nph
        sigline([nck(ii,1),nck(ii,2)],num2str(stats{ii}.p(2)),my(c)+my(c)*0.05*ii);
    end
end

% Extended Data Figure 10f - plot number of saccades
xvs = 0.35; %x-value scale that is used to space the data
pos = linspace(-xvs,xvs,length(filenames));
jitter = (pos(2)-pos(1))/2;
my = zeros(1,size(pTypes{1},1));
clear bl h hlegend
f1=figure('WindowState','maximized');hold on

i = find(contains(pTypes{1},'10NG'));
col = cReg(1,:);

for kk = 1 : length(genotypes)
    bl{kk} = SPKnum{1,kk}(:,i);
    maxbl = max(bl{kk});
    my(i) = max([my(i),maxbl]);

    if ~isempty(bl{kk})
        h{kk}=notBoxPlot(bl{kk},ones(1,length(bl{kk}))*...
            (i+pos(kk)),'jitter',jitter,'style','line');
        set(h{kk}.data,'color',[0 0 0],'MarkerFaceColor',...
            col,'MarkerSize',6,'LineWidth',0.5,...
            'Marker','o')
        set(h{kk}.mu,'color','k','MarkerFaceColor','k');
        set(h{kk}.sem,'visible','off');
        set(h{kk}.sd,'color','k');
    end
end

% Perform Mann-Whitney-Wilcoxon non parametric test and plot
% significance indicators if p<0.05
nck = nchoosek(1:kk,2);
for ii = 1 : size(nck,1)
    stats{ii} = mwwtest(bl{nck(ii,1)},bl{nck(ii,2)},0);
    if stats{ii}.p(2) < 0.05/nph
        sigline([i+pos(nck(ii,1)),i+pos(nck(ii,2))],num2str(stats{ii}.p(2)),my(i)+my(i)*0.05*ii);
    end
end
set(gca,'XTick',i+pos,'XTickLabel',gOrder,'Fontsize',18)
ylabel('Number of Saccades')
xlim([1-xvs*1.5 i+xvs*1.5])

%% Extended Data Figure 10e
% DNp15 hsFlp Kir
cond = {'Dark','Light'}; % experimental conditions
genotypes = {'Both','Left','None','Right'};
gTypes = {'Both','Left','None','Right'};
behDir = 'DNp15Kir\';
saveDir = [baseDir,'Analysis\DNp15Kir\']; % save directory
path = [baseDir,behDir,'Light'];
controls = [0 0 1 0];

[~,gi] = ismember(gOrder,gTypes); % plotting order indices


if exist([saveDir,'SaccadeVars.mat'],'file')
    load([saveDir,'SaccadeVars.mat'])
    [~,gi] = ismember(gOrder,gTypes); % plotting order indices
else
    cd(path);
    datedir = dir;
    filenames = {datedir(3:end).name};
    for k = 1 : length(filenames)
        pathF{k} = [path,'\',filenames{k},'\'];

        if ~isempty(pathF{k})
            [SPKVR{k}, ~, SPKVRDist{k},SPKVF{k}, pcav{k}, pTypes{k}, dt{k}] = GetSpikePCA(pathF{k}, params);
        end

        % visual stimuli order. Useful to plot in numerical order
        % (e.g. 1NG 5NG 10NG instead of 10NG 1NG 5NG)
        if strcmp('Dark',genotypes{k})
            seq{k} = 1;
        else
            seq{k} = [2,1];
        end
    end
end

% Plot saccade direction bias
fs = 14; % font size for plot labels
alpha = 0.7; % transparency of the dots
ms=32; % dot size of the mean
limsX = [0.3 length(genotypes)+0.7]; % axis limits
limsY = [-0.7 0.7];
% limsY = [-1 1];
scatterFactor = 0.16; % add scatter (fraction between 0 and 1) to the dots

switch cond{1}
    case 'Dark'
        xTicks = 1:length(genotypes);
        colors = [22,0,22]/256; conds = 1;
    case {'Light','Light Bilateral'}
        conds = seq{1};
        xTicks = 1:length(genotypes);
end

if contains(behDir,'DNp15')
    col = cReg(1,:);
elseif contains(behDir,'bIPS')
    col = cReg(2,:);
else
    col = [0,0,0];
end

my = zeros(1,length(conds));
SPKbinary = cell(size(SPKVR));
SPKnum = cell(size(SPKVR));
f=figure('WindowState','maximized');

c = find(contains(pTypes{1},'10NG'));
hold on;
for g = 1:length(genotypes) % for each genotype
    nFly = size(SPKVR{1,gi(g)},1);
    SPKbinary{1,g} = nan(nFly,length(cond(1)));
    for fly = 1:nFly % for each fly
        spikes = SPKVR{1,gi(g)}{fly,c};
        sbinary = nan(1,size(spikes,2));
        sbinary(max(spikes)<max(abs(spikes))) = -1;
        sbinary(max(spikes)==max(abs(spikes))) = 1;
        sbinary(spikes>0) = 1;
        sbinary(spikes<0) = -1;
        SPKbinary{1,g}(fly,c) = sum(sbinary)/length(sbinary);
        SPKnum{1,g}(fly,c) = length(sbinary);
    end
    vp = SPKbinary{1,g}(:,c); % value to plot
    vpn = SPKnum{1,g}(:,c); % num of events (for weighted vp)
    % remove flies without any saccades
    nd = sum(vpn==0); % number of discarded flies
    vp(vpn==0) = []; vpn(vpn==0) = [];

    Wmean = nansum((vp.*vpn)/sum(vpn));
    Wvar = nansum((((vp - Wmean).^2).*vpn)/sum(vpn));
    Wstd = sqrt(Wvar);

    sc = (rand(size(vp,1),1)-0.5)*scatterFactor;
    x = (ones(size(vp,1),1)*g) + sc;
    % plot the data, add mean and std error bars
    scatter(x,vp,vpn,...
        'MarkerFaceColor',col,'MarkerEdgeColor','none',...
        'MarkerFaceAlpha',alpha);

    plot(g,Wmean,'k.','MarkerSize',ms)
    eh=errorbar(g,Wmean,(1.960*Wstd/sqrt(length(vp))));
    set(eh,'LineWidth',2,'Color','k')

    maxbl = nanmax(SPKbinary{1,g}(:,c));
    my(c) = max(my(c),maxbl);
end
% set axis limits
xlim(limsX);
ylim(limsY);
ylabel('Saccade bias')
set(gca,'XTick',xTicks,'XTickLabel',gOrder,'FontSize',fs,...
    'XTickLabelRotation',45)
% add dashed line at 0
line(limsX,[0 0],'Color','k','LineWidth',0.5,'LineStyle','--')

nck = nchoosek(1:g,2);
nph = size(nck,1); % number of comparisons made (used for adjusting the
% minimum P value required for significance - a.k.a. Bonferroni correction)
for ii = 1 : size(nck,1)
    pair1 = SPKbinary{1,nck(ii,1)}(:,c);
    pair2 = SPKbinary{1,nck(ii,2)}(:,c);
    pair1(isnan(pair1)) = []; pair2(isnan(pair2)) = [];
    stats{ii} = mwwtest(pair1,pair2,0);
    if stats{ii}.p(2) < 0.05/nph
        sigline([nck(ii,1),nck(ii,2)],num2str(stats{ii}.p(2)),my(c)+my(c)*0.05*ii);
    end
end


%% Extended Data Figure 10g
% Folder with the preprocessed behavior data

% bIPS TrpA1
behDir = 'bIPSTrpA1\';
% experimental condition
cond = {'Light'};
genotypes = {'Empty Split Gal4 TrpA1','R25F07AD VT014724DBD TrpA1',...
    'VT014724AD VT039047DBD TrpA1','VT039047AD VT019845DBD TrpA1'};
gTypes = {'Empty','bIPS1','bIPS2','bIPS3'}; % short names for plot labels
saveDir = [baseDir,'Analysis\bIPSTrpA1\'];
path = [baseDir,behDir,'Light Bilateral'];

if exist([saveDir,'SaccadeVars.mat'],'file')
    load([saveDir,'SaccadeVars.mat'])
else
    cd(path);
    datedir = dir;
    filenames = {datedir(3:end).name};
    for k = 1 : length(filenames)
        pathF{k} = [path,'\',filenames{k},'\'];
        if ~isempty(pathF{k})
            [SPKVR{k}, ~, SPKVRDist{k},SPKVF{k}, pcav{k}, pTypes{k}, dt{k}] = GetSpikePCA(pathF{k}, params);
        end
        % visual stimuli order. Useful to plot in numerical order
        % (e.g. 1NG 5NG 10NG instead of 10NG 1NG 5NG)
        if strcmp('Dark',genotypes{k})
            seq{k} = 1;
        else
            seq{k} = [2,1];
        end
    end
end

xvs = 0.35; %x-value scale that is used to space the data
pos = linspace(-xvs,xvs,length(filenames));
jitter = (pos(2)-pos(1))/2;
my = zeros(1,size(pTypes{1},1));
clear bl h hlegend
f1=figure('WindowState','maximized');hold on

i = find(contains(pTypes{1},'10NG'));


for kk = 1 : length(genotypes)
    bl{kk} = SPKnum{1,kk}(:,i);
    maxbl = max(bl{kk});
    my(i) = max([my(i),maxbl]);

    if ~isempty(bl{kk})
        h{kk}=notBoxPlot(bl{kk},ones(1,length(bl{kk}))*...
            (i+pos(kk)),'jitter',jitter,'style','line');
        if kk==1
            col = [0,0,0];
        else
            col = cReg(2,:);
        end
        set(h{kk}.data,'color',[0 0 0],'MarkerFaceColor',...
            col,'MarkerSize',6,'LineWidth',0.5,...
            'Marker','o')
        set(h{kk}.mu,'color','k','MarkerFaceColor','k');
        set(h{kk}.sem,'visible','off');
        set(h{kk}.sd,'color','k');
    end
end

% Perform Mann-Whitney-Wilcoxon non parametric test and plot
% significance indicators if p<0.05
nck = nchoosek(1:kk,2);
for ii = 1 : size(nck,1)
    stats{ii} = mwwtest(bl{nck(ii,1)},bl{nck(ii,2)},0);
    if stats{ii}.p(2) < 0.05/nph
        sigline([i+pos(nck(ii,1)),i+pos(nck(ii,2))],num2str(stats{ii}.p(2)),my(i)+my(i)*0.05*ii);
    end
end
set(gca,'XTick',i+pos,'XTickLabel',gTypes,'Fontsize',18)
ylabel('Number of Saccades')
xlim([1-xvs*1.5 i+xvs*1.5])