% This script takes pre-processed imaging data from each fly|neuron|ROI,
% aligns them in a common template and plots the data shown in the paper.

clearvars
close all

%% Key Abbreviations
% -------------------------------------------------------------------------
% auc - Area under the curve
% CCW - Counter-clockwise
% CW - Clockwise
% DS - Direction selective
% DSI - Direction selectivity index
% DFF - Delta F over F, one of the primary measures of neural Ca2+ activity
% F - Raw fluorescence over time obtained from a given ROI
% ND - Non-preferred Direction (motion direction where the neuron minimally
% responds, is hyperpolarized or has a calcium response below its baseline)
% NPD - Same as ND, Non-preferred Direction
% OF - Optic flow
% PD - Preferred Direction (motion direction where the neuron is maximally
% depolarized)
% ROI - Region of Interest
% ST - Speed tuning (refers to protocols where a particular pattern of OF
% was presented to the fly at various speeds)
% YF - Yaw Flipped: Binocular Front-to-Back (bFtoB) and Back-to-Front
% (bBtoF) horizontal yaw motion
% -------------------------------------------------------------------------

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
dataDir = fullfile(repoDir, 'data', 'Preprocessed Data\');
saveDir = fullfile(repoDir, 'results', 'matlab_outputs\paper_figures\');

% Create the results directory if it doesn't exist
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

% Option flag to save figures (1 to save, 0 to just display)
SaveFigures = 1;

% Display the paths for verification
fprintf('Repository Directory: %s\n', repoDir);
fprintf('Data Directory: %s\n', dataDir);
fprintf('Results Directory: %s\n', saveDir);

% remove unnecessary variables
clearvars homeDir mainDir

%% Load the data

% find out how many .mat files exists in the analysis directory
numFiles=length(subdir(fullfile(dataDir,'*.mat')));
progressbar('Loading all sessions ...') % initialize progress bar

% All data will be stored in a structure file, whose size will depend on
% the number of groups and number of imaging sessions within each group
data = struct([]);

% Parse all folders in the data directory. It has the following format:
% dataDir > neurons > genotypes > date of imaging > fly number > .mat files
neurons = dir(dataDir); neurons(1:2)=[];
neuronNames = {neurons.name};
for n = 1:length(neurons) % for each neuron
    gDir = [dataDir,neurons(n).name];
    genotypes = dir(gDir);
    % remove the first 2 entries (windows has 2 extra dirs in each folder)
    genotypes(1:2)=[]; % this is done for every subfolder below
    if ~isempty(genotypes)
        for g = 1:length(genotypes) % for each genotype within a neuron
            dDir = [gDir,'\',genotypes(g).name];
            dates = dir(dDir); dates(1:2) = [];
            for d = 1:length(dates) % for each day within a genotype
                fDir = [dDir,'\',dates(d).name];
                flies = dir(fDir); flies(1:2) = [];
                for f = 1:length(flies) % for each fly within a day
                    sDir = [fDir,'\',flies(f).name,'\'];
                    sessions = dir(sDir); sessions(1:2) = [];
                    for s = 1:length(sessions) % for each session per fly
                        fNum = size(data,2)+1; % iterate for each .mat file
                        progressbar(fNum/numFiles)% update the progress bar

                        % load .mat file to a temporary folder and only
                        % pass the relevant fields to the data structure.
                        % If you are processing many flies and pass
                        % everything directly, you may have RAM overflow
                        tempData = load([sDir,sessions(s).name]);
                        if ishandle(2)
                            disp([sDir,sessions(s).name])
                            close(2)
                        end
                        % pass speed tuning data as well, if it exists
                        if isfield(tempData,'DFF_ST')
                            % extract motion direction information
                            fnc = tempData.HeaderFile.AcqInfo...
                                .FunctionNameMat{1}{1};
                            if contains(fnc,'-CCW-')
                                direction = 'CCW';
                            else
                                if contains(fnc,'-CW-')
                                    direction = 'CW';
                                else
                                    error(['Speed tuning function',...
                                        'direction cannot be detected!'])
                                end
                            end
                            % append motion direction to the pattern names
                            for po = 1:length(tempData.PatsOrdered)
                                tempStr = tempData.PatsOrdered{po};
                                usi = strfind(tempStr,'_');
                                tempStr = [tempStr(1:usi(end)),...
                                    direction,tempStr(usi(end):end)];
                                tempData.PatsOrdered{po} = tempStr;
                            end
                            % save the relevant variables to data structure
                            data(fNum).data = struct(...
                                'F_OF',tempData.F_All,...
                                'F_ST',tempData.F_ST,...
                                'DFF_OF',tempData.DFF_All,...
                                'DFF_ST',tempData.DFF_ST,...
                                'ROIsToBatchAnalyze',...
                                tempData.ROIsToBatchAnalyze,...
                                'PatsOrdered',{tempData.PatsOrdered});
                        else
                            data(fNum).data = struct(...
                                'F_OF',tempData.F_All,...
                                'DFF_OF',tempData.DFF_All,...
                                'ROIsToBatchAnalyze',...
                                tempData.ROIsToBatchAnalyze,...
                                'PatsOrdered',{tempData.PatsOrdered});
                        end
                        % pass the relevant meta-data
                        data(fNum).genotype = genotypes(g).name;
                        if contains(tempData.neuron,'H2in')
                            tempData.neuron = ...
                                strrep(tempData.neuron,'H2in','H2rn');
                        end
                        data(fNum).neuron = tempData.neuron;
                        data(fNum).date = dates(d).name;
                        data(fNum).hemisphere = tempData.hemi;
                        data(fNum).region = tempData.neurite;
                        data(fNum).somaSide = tempData.soma;
                        data(fNum).prepType = tempData.prepType;

                        % unique identifier for each neurite. Format is
                        % Date-FlyNumber-SomaLocation-Neuron-Neurite
                        % this "barcode" will allow to pool all sessions
                        % coming from the same fly>neuron>neurite
                        data(fNum).barcode = [dates(d).name,...
                            flies(f).name,data(fNum).somaSide,...
                            data(fNum).neuron,data(fNum).region];

                        % turning frame rates into integers for resample
                        % function. NOTE: we are not multiplying by 10^5
                        % because then the integers become too big for the
                        % resample function. instead, we settle for 5-digit
                        % frame rates. Resampling is done later below
                        data(fNum).frameRate = ...
                            round(tempData.FrameRate*1000);
                    end
                end
            end
        end
    end
end
% remove unnecessary variables
clearvars d dDir dates f fDir flies fNum g gDir genotypes n neurons s ...
    sessions tempData sDir fnc direction tempStr usi po

%% Pool imaging data from multiple sessions done on the same neurite

% The term neurite refers to a particular branch (axon/dendrite/mixed) of a
% single cell. In some flies multiple cell types were imaged and in some
% flies axons and dendrites of the same cell were imaged in separate
% sessions. Therefore the best way to organize (and pool) data is done on a
% neurite by neurite basis. That's why the unique "barcode" per fly
% includes the cell type and the imaged neurite.
%
% Many ROIs were imaged for multiple protocols (rotational vs translational
% OF sensitivities, speed tuning, roll axis sensitivity etc). Here we pool
% all of that data into two variables per neurite.
%
% Speed tuning protocols are shorter (6 seconds) than the rotational and
% translational OF panel stimuli (13 seconds), which makes it impractical
% to pool and analyze together, even if they share some of the same pattern
% and speed combinations. Therefore speed tuning protocols are separated
% from the main OF.data.dFF matrix and saved in the OF.data.dFFST variable.
% Stimulus names are also separated accordingly, within OF.data variable.

OF=data([]);

% Each protocol may have different sets of visual stimuli shown. These
% variables will serve as a template to visualize all of the stimuli shown,
% and to sort imaging data from each fly to this template later on

Stimuli = {}; % stimuli used in OF sensitivity and receptive field mapping
StimuliST = {}; % stimuli used in speed tuning

neurites = unique({data.barcode});

% The lowest imaging frame rate used in all of the experiments
minFrameRate = min([data.frameRate]);

% initialize progress bar and the progress counter
prog = 0; progressbar('Pooling all sessions from each neurite...')

for n = 1:length(neurites) % for each neurite imaged

    % check if the same ROI has been used in other protocols
    ni=find(strcmp(neurites{n},{data.barcode})); %neurite index. it will
    % include all the indices (row number in data) where this neurite was
    % imaged

    fNum = size(OF,2)+1; % iterate for each new neurite

    for d = 1:size(ni,2) % for every protocol where this neurite was imaged
        prog=prog+1; progressbar(prog/numFiles) % update the progress bar
        if d==1 % initialize OF.data structure
            OF(fNum) = data(ni(d)); OF(fNum).data = struct('F',[],...
                'DFF',[],'FST',[],'DFFST',[],...
                'Stimuli',[],'StimuliST',[],...
                'F_All_concatenated',[],'DFF_All_concatenated',[]);
        end
        if isfield(data(ni(d)).data,'DFF_OF') % DFF_All (renamed to DFF_OF
            % in this script during data loading phase) inside each file is
            % a 4D matrix where deltaF/F calcium responses are stored in
            % the following format: (FrameNum,TrialNum,VisualStimulus,ROIs)

            % subset of ROIs to be analyzed
            r = data(ni(d)).data.ROIsToBatchAnalyze;

            if isfield(data(ni(d)).data,'DFF_ST') % this variable should
                % only exist (and saved) in flies where we performed speed
                % tuning experiments

                % Downsample data to lowest existing frame rate. This is
                % done to be able to align time series data coming from
                % experiments with different acquisition frame rates
                ResampledData = resampleImagingData(...
                    data(ni(d)).data.DFF_ST(:,:,:,r),...
                    data(ni(d)).frameRate,minFrameRate); % DFF
                ResampledDataF = resampleImagingData(...
                    data(ni(d)).data.F_ST(:,:,:,r),...
                    data(ni(d)).frameRate,minFrameRate); % F

                if isempty(OF(fNum).data.DFFST) % add the first data point
                    OF(fNum).data.DFFST = ResampledData;
                    OF(fNum).data.FST = ResampledDataF;
                    OF(fNum).data.StimuliST = data(ni(d)).data.PatsOrdered;
                else
                    % if OF.data.DFFST already has some speed tuning data,
                    % detect dimensions of OF.data with the incoming DFF
                    % matrix and concatenate them together, then update the
                    % OF.data.DFFST variable
                    dimOFST = size(OF(fNum).data.DFFST);
                    dimdataST = size(ResampledData);
                    % if there's only one ROI, add 4th dimension explicitly
                    if length(dimOFST)==3
                        dimOFST = [dimOFST,1]; %#ok<*AGROW>
                        % comment above is used to suppress warnings
                    end
                    if length(dimdataST)==3
                        dimdataST = [dimdataST,1];
                    end
                    dimnewST = max([dimOFST;dimdataST]);
                    dimnewST(3) = dimOFST(3)+dimdataST(3);
                    newMatST = nan(dimnewST);
                    newMatST(1:dimOFST(1),1:dimOFST(2),1:dimOFST(3),...
                        1:dimOFST(4)) = OF(fNum).data.DFFST;
                    newMatST(1:dimdataST(1),1:dimdataST(2),...
                        (1+dimOFST(3)):end,1:dimOFST(4)) = ResampledData;
                    OF(fNum).data.DFFST = newMatST;

                    % same concatenation for the F time series
                    dimOFSTF = size(OF(fNum).data.FST);
                    dimdataSTF = size(ResampledDataF);
                    if length(dimOFSTF)==3
                        dimOFSTF = [dimOFSTF,1]; %#ok<*AGROW>
                        % comment above is used to suppress warnings
                    end
                    if length(dimdataSTF)==3
                        dimdataSTF = [dimdataSTF,1];
                    end
                    dimnewSTF = max([dimOFSTF;dimdataSTF]);
                    dimnewSTF(3) = dimOFSTF(3)+dimdataSTF(3);
                    newMatSTF = nan(dimnewSTF);
                    newMatSTF(1:dimOFSTF(1),1:dimOFSTF(2),1:dimOFSTF(3),...
                        1:dimOFSTF(4)) = OF(fNum).data.FST;
                    newMatSTF(1:dimdataSTF(1),1:dimdataSTF(2),...
                        (1+dimOFSTF(3)):end,1:dimOFSTF(4)) = ResampledData;
                    OF(fNum).data.FST = newMatSTF;

                    % stimulus info
                    OF(fNum).data.StimuliST = [OF(fNum).data.StimuliST,...
                        data(ni(d)).data.PatsOrdered];
                end
            else
                % downsample data to lowest existing frame rate
                ResampledData = resampleImagingData(...
                    data(ni(d)).data.DFF_OF(:,:,:,r),...
                    data(ni(d)).frameRate,minFrameRate); %DFF
                ResampledDataF = resampleImagingData(...
                    data(ni(d)).data.F_OF(:,:,:,r),...
                    data(ni(d)).frameRate,minFrameRate); %F
                % Ignore if the protocol is roll only (roll OF tuning)
                if ~all(contains(data(ni(d)).data.PatsOrdered,'DotRollInv'))
                    if isempty(OF(fNum).data.DFF) % add first data point
                        OF(fNum).data.DFF = ResampledData;
                        OF(fNum).data.F = ResampledDataF;
                        % Little hack to give a unique name to the latest
                        % flicker control in the trans+rot compound stimuli
                        % Otherwise it interferes with sorting the
                        % similarly named Trans2-14Rot1 stimuli
                        data(ni(d)).data.PatsOrdered = ...
                            strrep(data(ni(d)).data.PatsOrdered,...
                            'Flicker_Trans2-14Rot1','FlickerTransRot');
                        OF(fNum).data.Stimuli = ...
                            data(ni(d)).data.PatsOrdered;
                    else
                        % if OF.data.DFF already has some data, detect the
                        % dimensions of OF.data.DFF and the incoming DFF
                        % matrix, concatenate them together, then update
                        % the OF.data.DFF variable
                        dimOF = size(OF(fNum).data.DFF);
                        dimdata = size(ResampledData);
                        % if there's one ROI, add 4th dimension explicitly
                        if length(dimOF)==3
                            dimOF = [dimOF,1];
                        end
                        if length(dimdata)==3
                            dimdata = [dimdata,1];
                        end
                        dimnew = max([dimOF;dimdata]);
                        dimnew(3) = dimOF(3)+dimdata(3);
                        tempMat = nan(dimnew);
                        tempMat(1:dimOF(1),1:dimOF(2),1:dimOF(3),...
                            1:dimOF(4))=OF(fNum).data.DFF;
                        tempMat(1:dimdata(1),1:dimdata(2),...
                            (1+dimOF(3)):end,1:dimOF(4)) = ResampledData;
                        OF(fNum).data.DFF = tempMat;

                        % same concatenation for the F time series
                        dimOFF = size(OF(fNum).data.F);
                        dimdataF = size(ResampledDataF);
                        if length(dimOFF)==3
                            dimOFF = [dimOFF,1];
                        end
                        if length(dimdataF)==3
                            dimdataF = [dimdataF,1];
                        end
                        dimnewF = max([dimOFF;dimdataF]);
                        dimnewF(3) = dimOFF(3)+dimdataF(3);
                        tempMatF = nan(dimnewF);
                        tempMatF(1:dimOFF(1),1:dimOFF(2),...
                            1:dimOFF(3),1:dimOFF(4)) = OF(fNum).data.F;
                        tempMatF(1:dimdataF(1),1:dimdataF(2),...
                            (1+dimOFF(3)):end,1:dimOFF(4))=ResampledDataF;
                        OF(fNum).data.F = tempMatF;


                        % Little hack to give a unique name to the latest
                        % flicker control in the trans+rot compound stimuli
                        % Otherwise it interferes with sorting the
                        % similarly named Trans2-14Rot1 stimuli
                        data(ni(d)).data.PatsOrdered = ...
                            strrep(data(ni(d)).data.PatsOrdered,...
                            'Flicker_Trans2-14Rot1','FlickerTransRot');
                        % stimulus info
                        OF(fNum).data.Stimuli = [OF(fNum).data.Stimuli,...
                            data(ni(d)).data.PatsOrdered];

                        % In some cases, some of the identical stimuli have
                        % been used across multiple protocols (e.g.
                        % bilateral yaw in rotational OF protocol and HS
                        % receptive field mapping) The code below detects
                        % such duplicates in the final OF.data matrix and
                        % pool identical stimuli together into the same
                        % columns in OF.data.DFF and OF.data.F variables
                        [~,fid,uid]=unique(OF(fNum).data.Stimuli,'stable');
                        if length(uid)>max(uid) % if there are duplicates
                            nuid = []; % initialize non-unique indices var.
                            % get the current dimensions of DFF matrix
                            dimOF = size(OF(fNum).data.DFF);
                            % if there's one ROI, add 4th dim. explicitly
                            if length(dimOF)==3
                                dimOF = [dimOF,1];
                            end
                            % initialize temp. matrices without duplicates
                            tempMat = nan(dimOF);
                            tempMatF = nan(dimOF);
                            for i = 1:length(uid) % for each index of new
                                % matrices
                                if ismember(i,fid) % if non-duplicate (1st)
                                    dimTemp = size(OF(fNum).data.DFF);
                                    % if there's one ROI, add 4th dim.
                                    if length(dimTemp)==3
                                        dimTemp = [dimTemp,1];
                                    end
                                    tempMat(1:dimTemp(1),...
                                        1:dimTemp(2),...
                                        i,1:dimTemp(4)) = ...
                                        OF(fNum).data.DFF(:,:,i,:);
                                    tempMatF(1:dimTemp(1),...
                                        1:dimTemp(2),...
                                        i,1:dimTemp(4)) = ...
                                        OF(fNum).data.F(:,:,i,:);
                                else % if duplicate stimulus (2nd or more)
                                    % calculate new dimentions and update
                                    % the temporary matrices
                                    dim2 = ...
                                        size(OF(fNum).data.DFF(:,:,i,:),2);
                                    dimOFnew = dimOF;
                                    dimOFnew(2) = dimOF(2)+dim2;
                                    tempMat2 = nan(dimOFnew);
                                    tempMat2(1:dimOF(1),...
                                        1:dimOF(2),...
                                        1:dimOF(3),...
                                        1:dimOF(4)) = tempMat;
                                    tempMat2(1:dimOFnew(1),...
                                        dimOF(2)+1:dimOFnew(2),...
                                        uid(i),1:dimOFnew(4)) = ...
                                        OF(fNum).data.DFF(:,:,i,:);
                                    tempMat = tempMat2;
                                    tempMatF2 = nan(dimOFnew);
                                    tempMatF2(1:dimOF(1),...
                                        1:dimOF(2),...
                                        1:dimOF(3),...
                                        1:dimOF(4)) = tempMatF;
                                    tempMatF2(1:dimOFnew(1),...
                                        dimOF(2)+1:dimOFnew(2),...
                                        uid(i),1:dimOFnew(4)) = ...
                                        OF(fNum).data.F(:,:,i,:);
                                    tempMatF = tempMatF2;
                                    dimOF = dimOFnew;
                                    nuid = [nuid,i]; % non-unique indices
                                end
                            end
                            % once you concatenate all of the duplicate
                            % columns, the temporary matrices will be full
                            % of unnecessary NaN columns. We remove those
                            % NaNs below, before passing the data onto the
                            % OF.data variable.

                            % Remove duplicate column(s) from temp matrix
                            tempMat(:,:,nuid,:) = [];
                            tempMatF(:,:,nuid,:) = [];
                            % remove NaN trials from each stimulus
                            for i = 1:size(tempMat,3) % for each stimulus
                                % get the DFF data
                                tempMat2 = squeeze(tempMat(:,:,i));
                                tempMatF2 = squeeze(tempMatF(:,:,i));
                                % remove NaN columns
                                tempMat3 = rmmissing(tempMat2,2);
                                tempMatF3 = rmmissing(tempMatF2,2);
                                if i==1 % if it is the first entry
                                    % generate temporary matrix with the
                                    % new dimensions
                                    dimOFnew = dimOF;
                                    dimOFnew(2) = size(tempMat3,2);
                                    tempMat4 = nan(dimOFnew);
                                    tempMat4(:,:,i,:) = tempMat3;
                                    tempMatF4 = nan(dimOFnew);
                                    tempMatF4(:,:,i,:) = tempMatF3;
                                else
                                    % num of trials per subsequent stimulus
                                    dimTrials = size(tempMat3,2);
                                    if dimOFnew(2)>= dimTrials % if smaller
                                        % or equal to new matrix dimensions
                                        % simply add them to the new matrix
                                        tempMat4(:,1:dimTrials,i,:) = ...
                                            tempMat3;
                                        tempMatF4(:,1:dimTrials,i,:) = ...
                                            tempMatF3;
                                    else % if there are more trials than
                                        % current matrix dimensions, update
                                        % the temporary matrix with the new
                                        % dimensions
                                        dimOFnew2 = dimOFnew;
                                        % total num of trials = new maximum
                                        dimOFnew2(2) = size(tempMat3,2);
                                        % tempMat5 to hold new dimentions
                                        tempMat5 = nan(dimOFnew2);
                                        tempMatF5 = nan(dimOFnew2);
                                        % add the previous data to temp5
                                        tempMat5(1:dimOFnew(1),...
                                            1:dimOFnew(2),1:dimOFnew(3),...
                                            1:dimOFnew(4)) = tempMat4;
                                        tempMatF5(1:dimOFnew(1),...
                                            1:dimOFnew(2),1:dimOFnew(3),...
                                            1:dimOFnew(4)) = tempMatF4;
                                        % add the newest stimulus data
                                        tempMat5(:,:,i,:) = tempMat3;
                                        tempMatF5(:,:,i,:) = tempMatF3;
                                        % update temp mat
                                        tempMat4 = tempMat5;
                                        tempMatF4 = tempMatF5;
                                        % update dimensions
                                        dimOFnew = dimOFnew2;
                                    end
                                end
                            end
                            % update OF.data
                            OF(fNum).data.DFF = tempMat4;
                            OF(fNum).data.F = tempMatF4;
                            % remove duplicates from the Stimuli name array
                            OF(fNum).data.Stimuli(nuid) = [];
                        end
                    end
                end
            end
        else
            warning(['no data gathered from ',data(ni(d)).barcode,'-',...
                num2str(d),', check associated .mat file'])
        end
    end

    % Remove the zeroes that were added during matrix concatenations
    OF(fNum).data.DFF(OF(fNum).data.DFF==0) = nan;
    OF(fNum).data.DFFST(OF(fNum).data.DFFST==0) = nan;
    OF(fNum).data.F(OF(fNum).data.F==0) = nan;
    OF(fNum).data.FST(OF(fNum).data.FST==0) = nan;

    % Stimuli and StimuliST will contain every protocol used in the entire
    % collection of flies, such that each fly/neurite can be mapped onto a
    % subset of the Stimuli and StimuliST cell arrays later on.
    % Append new visual stimuli that don't exist in Stimuli & StimuliST
    if ~isempty(OF(fNum).data.DFF) % If neurite has imaging data
        if isempty(Stimuli)
            Stimuli = OF(fNum).data.Stimuli;
        else
            for s = 1:length(OF(fNum).data.Stimuli)
                if ~sum(strcmp(Stimuli,OF(fNum).data.Stimuli{s}))
                    tempStr = [Stimuli OF(fNum).data.Stimuli{s}];
                    Stimuli = tempStr;
                end
            end
        end
    else
        warning([OF(fNum).barcode,...
            ' does not have imaging data, check .mat file!'])
    end
    if ~isempty(OF(fNum).data.DFFST) % If neurite has speed tuning data
        if isempty(StimuliST)
            StimuliST = OF(fNum).data.StimuliST;
        else
            for s = 1:length(OF(fNum).data.StimuliST)
                if ~sum(strcmp(StimuliST,OF(fNum).data.StimuliST{s}))
                    tempStr = [StimuliST OF(fNum).data.StimuliST{s}];
                    StimuliST = tempStr;
                end
            end
        end
    end

    % find number of ROIs per neurite
    numROIs = max([size(OF(n).data.DFF,4) size(OF(n).data.DFFST,4)]);
    % concatenate all imaging data per ROI into a single vector
    for ROI = 1:numROIs
        if ~isempty(OF(fNum).data.DFF)
            tempdFFPerROI = squeeze(OF(n).data.DFF(:,:,:,ROI));
            tempdFFPerROI = ...
                reshape(tempdFFPerROI,[1,numel(tempdFFPerROI)]);
            tempdFFPerROI(isnan(tempdFFPerROI)) = [];

            tempFPerROI = squeeze(OF(n).data.F(:,:,:,ROI));
            tempFPerROI = reshape(tempFPerROI,[1,numel(tempFPerROI)]);
            tempFPerROI(isnan(tempFPerROI)) = [];
        else
            tempdFFPerROI = [];
            tempFPerROI = [];
        end
        if ~isempty(OF(fNum).data.DFFST)
            tempdFFSTPerROI = squeeze(OF(n).data.DFFST(:,:,:,ROI));
            tempdFFSTPerROI = ...
                reshape(tempdFFSTPerROI,[1,numel(tempdFFSTPerROI)]);
            tempdFFSTPerROI(isnan(tempdFFSTPerROI)) = [];

            tempFSTPerROI = squeeze(OF(n).data.FST(:,:,:,ROI));
            tempFSTPerROI = ...
                reshape(tempFSTPerROI,[1,numel(tempFSTPerROI)]);
            tempFSTPerROI(isnan(tempFSTPerROI)) = [];
        else
            tempdFFSTPerROI = [];
            tempFSTPerROI = [];
        end
        OF(n).data.DFF_All_concatenated(:,:,ROI) = [tempdFFPerROI,...
            tempdFFSTPerROI];
        OF(n).data.F_All_concatenated(:,:,ROI) = [tempFPerROI,...
            tempFSTPerROI];
    end
end

% Order the Stimuli names alphabetically
Stimuli = sort(Stimuli);

% StimuliST should already be sorted numerically and shouldn't be resorted

% Transpose the cell arrays for easier plotting of stimuli names in the
% command window
Stimuli = Stimuli';
StimuliST = StimuliST';

% remove unnecessary variables
clearvars d dimdata dimnew dimOF fNum groupID n ni nu nuIdx r dashes ...
    ResampledData s tempMat tempStr neurites prog numFiles dimdataST ...
    dimnewST dimOFST dimtrials ROI newMatST newOF nuOF I numbers s i ...
    tempdFFPerROI tempdFFSTPerROI numROIs numFiles dimOFnew dimOFnew2 ...
    ResampledDataF dimOFSTF dimdataSTF dimnewSTF newMatSTF dimOFF dim2 ...
    dimdataF dimnewF tempMatF newOFF nuOFF dimtrialsF tempFPerROI fid ...
    uid nuid IC tempFSTPerROI tempMatF2 tempMatF3 tempMatF4 tempMatF5 ...
    tempMat2 tempMat3 tempMat4 tempMat5 dimTrials dimTemp

%% Sort imaging data with the same stimulus order for all flies

% Here we reorder all of the imaging data within the variables OF.data.F,
% OF.data.DFF and OF.data.DFFST for each neurite onto a master template.
% Stimulus order is stored within the Stimuli and StimuliST variables.

for n = 1:length(OF) % for each neurite
    % create a temporary array with columns for all the stimuli presented
    tempdff = nan(size(OF(n).data.DFF,1),size(OF(n).data.DFF,2),...
        length(Stimuli),size(OF(n).data.DFF,4));
    tempf = nan(size(OF(n).data.F,1),size(OF(n).data.F,2),...
        length(Stimuli),size(OF(n).data.F,4));

    % allocate the data from DFF and F arrays into the correct columns of
    % the tempdff and tempf temporary arrays
    for s = 1:length(OF(n).data.Stimuli) % for each visual stimulus
        tempdff(:,:,contains(Stimuli,OF(n).data.Stimuli{s}),:) = ...
            OF(n).data.DFF(:,:,s,:);
        tempf(:,:,contains(Stimuli,OF(n).data.Stimuli{s}),:) = ...
            OF(n).data.F(:,:,s,:);
    end

    % rewrite the imaging data
    OF(n).data.DFF = tempdff;
    OF(n).data.F = tempf;

    % remove the outdated Stimuli field from structure
    OF(n).data=rmfield(OF(n).data,'Stimuli');

    % do the same reordering with the speed tuning data, if it exists
    if ~isempty(OF(n).data.DFFST)
        tempdff = nan(size(OF(n).data.DFFST,1),size(OF(n).data.DFFST,2),...
            length(StimuliST),size(OF(n).data.DFFST,4));
        tempf = nan(size(OF(n).data.FST,1),size(OF(n).data.FST,2),...
            length(StimuliST),size(OF(n).data.FST,4));
        for s = 1:length(OF(n).data.StimuliST) % for each visual stimulus
            tempdff(:,:,contains(StimuliST,OF(n).data.StimuliST{s}),:)=...
                OF(n).data.DFFST(:,:,s,:);
            tempf(:,:,contains(StimuliST,OF(n).data.StimuliST{s}),:)=...
                OF(n).data.FST(:,:,s,:);
        end
        OF(n).data.DFFST = tempdff;
        OF(n).data.FST = tempf;
    end
    OF(n).data=rmfield(OF(n).data,'StimuliST');
end

clearvars n tempdff s tempf

%% Generate ipsi-, contra- and bi-lateral response matrices

% As of 210414 all OF stimuli contains the structure of 2-2-5-2-2 seconds
% (still - motion in one direction - still - reverse motion - still)
% except for speed tuning stimuli, which has motion only in one direction
% 2-2-2 sec (still - motion - still). Therefore, it is required to divide
% the stimuli into its separate motion "epochs" and allocate them into
% matrices that depict the soma location of recorded neurons (Left/Right).
% 240105: Translational + Rotational compound OF stimuli also follow the
% 2-2-2 sec structure, even though they are not speed tuning stimuli. We
% treat them slightly differently below.

% OF(:).data.DFF_1st and OF(:).data.DFF_2nd will contain DFF traces where
% the recordings from Left and Right cells are pooled together by flipping
% the Right cell recordings; OF(:).data.DFF_1st_Raw and
% OF(:).data.DFF_2nd_Raw will contain DFF traces that are not pooled for
% moments when you need to analyze Left and Right recordings separately.
% See comments later in this section for exactly how the recordings from
% Right cells are flipped.
% OF(:).data.DFF_1st_R and OF(:).data.DFF_2nd_R will contain pooled
% recordings, but for the Right hemisphere cells (used mainly for H2)
% OF(:).data.F_ variables are the same as DFF but for raw Fluorescence

% Separate visual stimuli indices based on stimulus laterality (stim shown
% on the Left/Right half of the LED arena or binocularly aka on both sides)
LeftIndices = contains(Stimuli,'_Left_');
RightIndices = contains(Stimuli,'_Right_');
BinocularIndices = ~LeftIndices & ~RightIndices;

StimuliLeft = Stimuli(LeftIndices);
StimuliRight = Stimuli(RightIndices);
StimuliBinocular = Stimuli(BinocularIndices);

% As of 221104 the default is to organize all the recordings and OF Stimuli
% presentation to match the responses recorded from Left hemisphere neurons
% Right hemisphere recordings therefore are flipped accordingly. The
% StimuliFlipped array instead organizes OF Stimuli for Right hemisphere
% neurons. EpochsFlipped indicates if the 1st Epoch remains the same (1) or
% is flipped (2) after transposing stimulus patterns to Right hemisphere
StimuliFlipped = cell(size(Stimuli));
EpochsFlipped = zeros(size(Stimuli));

% separate visual stimuli indices based on midline symmetry. As of 210414
% the visual stimuli that are symmetric at the visual display midline are:
% Pitch, Progressive, Lift and YawFlipped. Only Pitch and Progressive have
% trials where the pattern is shown only on one side of the LED arena

% indices of monocular presentation of midline symmetric OF stimuli
iSymUni = contains(StimuliLeft,{'Progressive','Pitch'});
% indices of binocular presentation of midline symmetric OF stimuli
iSymBilat = contains(StimuliBinocular,{'Progressive','Pitch'});
% indices of monocular presentation of midline asymmetric OF stimuli
iYawUni = contains(StimuliLeft,{'Yaw_'}); % Yaw index
iRollUni = contains(StimuliLeft,{'RollInv_'}); % Roll index
% indices of binocular presentation of midline asymmetric OF stimuli
iYawBilat = contains(StimuliBinocular,{'Yaw_'}); % Yaw index
iRollBilat = contains(StimuliBinocular,{'RollInv_'}); % Roll index
iSsBilat = contains(StimuliBinocular,{'Sideslip_'}); % Sideslip index
iYF = contains(StimuliBinocular,{'YawRightFlipped'}); % YawFlipped index
iLiftBilat = contains(StimuliBinocular,{'Lift'}); % Lift index

% Indices of OF patterns that are symmetric at the visual display midline
% (both monocular and binocular presentations are included)
iSym = contains(Stimuli,{'Progressive','Pitch','Lift','YawRightFlipped'});
% Indices of stimuli that are asymmetric at the visual display midline
iAsym = contains(Stimuli,{'Yaw_','Sideslip','RollInv'});
iProgSide45 = contains(Stimuli,{'ProgSide45'});
iProgSide315 = contains(Stimuli,{'ProgSide315'});
iTrans = contains(Stimuli,{'Trans'});
iTransCCW = contains(Stimuli,{'Trans'}) & contains(Stimuli,{'_CCW_'});
iTransCW = contains(Stimuli,{'Trans'}) & contains(Stimuli,{'_CW_'});

% Divide each stimulus between its two epochs, based on the fixed timing
% structure of the OF stimuli. (see description at the beginning of this
% section) As of 210414 --> 1st epoch: 0-6 seconds, 2nd epoch: 7-13 seconds
TimeStructure = cumsum([2 2 5 2 2]);
AddTime = 2; % Additional time before and after OF motion that is used to
% defining motion epochs and plotting calcium responses to a given visual
% stimulus (in seconds)

frameNum = size(OF(1).data.DFF,1); % number of frames per trial
% Define x-axis (Time) based on the final frame number after resampling
TimeAxis = linspace(0,frameNum/(minFrameRate/1000),frameNum);

% Define the time indices for baseline and motion epochs
blTime = [1 2]; % Time interval where the baseline mu and sigma are
% calculated for z-scoring. As of 210414 this is done between seconds [1-2]
% which is the time interval between one second after displaying the static
% pattern and the onset of the pattern motion
t=blTime(1); % Baseline start (in seconds)
[~,~,idx]=unique(abs(TimeAxis-t));
i_t0_start = find(idx==1); % Baseline start (in frames)
t=blTime(2); % Baseline end (in seconds)
[~,~,idx]=unique(abs(TimeAxis-t));
i_t0_end = find(idx==1); % Baseline end (in frames)
% start index for 1st epoch
t=TimeStructure(1)-AddTime;
[~,~,idx]=unique(abs(TimeAxis-t));
i_t1_start = find(idx==1);
% end index for 1st epoch
t=TimeStructure(2)+AddTime;
[~,~,idx]=unique(abs(TimeAxis-t));
i_t1_end = find(idx==1);
% start index for 2nd epoch
t=TimeStructure(3)-AddTime;
[~,~,idx]=unique(abs(TimeAxis-t));
i_t2_start = find(idx==1);
% end index for 2nd epoch
t=TimeStructure(4)+AddTime;
[~,~,idx]=unique(abs(TimeAxis-t));
i_t2_end = find(idx==1);
% All frames of baseline, the 1st and the 2nd epochs
t0 = i_t0_start:i_t0_end; % Baseline indices (in frames)
t1 = i_t1_start:i_t1_end; % Epoch 1 (in frames)
t2 = i_t2_start:i_t2_end; % Epoch 2 (in frames)
% If there are unintended extra frames, ignore them
if length(t1) > length(t2)
    t1 = t1(1:length(t2));
else
    if length(t1) < length(t2)
        t2 = t2(end-length(t1):end);
    end
end

% populate F and DFF matrices with baseline and movement epoch time series
for n = 1:length(OF) % for each neurite
    if ~isempty(OF(n).data.DFF) % if there is DFF data
        % if there is onlz shorter protocols, append the data matrix w/ NaN
        if size(OF(n).data.F,1) <= length(t1)
            OF(n).data.F((size(OF(n).data.F,1)+1):t2(end),:,:) = NaN;
            OF(n).data.DFF((size(OF(n).data.DFF,1)+1):t2(end),:,:) = NaN;
        end
        % initilaize F and DFF arrays for baseline and each movement epoch
        OF(n).data.F_1st = nan(length(t1),size(OF(n).data.F,2),...
            size(OF(n).data.F,3),size(OF(n).data.F,4));
        OF(n).data.F_2nd = OF(n).data.F_1st;
        OF(n).data.F_BL = nan(length(t0),size(OF(n).data.F,2),...
            size(OF(n).data.F,3),size(OF(n).data.F,4));

        OF(n).data.DFF_1st = nan(length(t1),size(OF(n).data.DFF,2),...
            size(OF(n).data.DFF,3),size(OF(n).data.DFF,4));
        OF(n).data.DFF_2nd = OF(n).data.DFF_1st;
        OF(n).data.DFF_BL = nan(length(t0),size(OF(n).data.DFF,2),...
            size(OF(n).data.DFF,3),size(OF(n).data.DFF,4));

        OF(n).data.F_1st_R = nan(length(t1),size(OF(n).data.F,2),...
            size(OF(n).data.F,3),size(OF(n).data.F,4));
        OF(n).data.F_2nd_R = OF(n).data.F_1st_R;
        OF(n).data.F_BL_R = nan(length(t0),size(OF(n).data.F,2),...
            size(OF(n).data.F,3),size(OF(n).data.F,4));

        OF(n).data.DFF_1st_R = nan(length(t1),size(OF(n).data.DFF,2),...
            size(OF(n).data.DFF,3),size(OF(n).data.DFF,4));
        OF(n).data.DFF_2nd_R = OF(n).data.DFF_1st_R;
        OF(n).data.DFF_BL_R = nan(length(t0),size(OF(n).data.DFF,2),...
            size(OF(n).data.DFF,3),size(OF(n).data.DFF,4));

        % Pass the data of midline symmetric OF stimuli
        OF(n).data.F_1st(:,:,iSym,:) = OF(n).data.F(t1,:,iSym,:);
        OF(n).data.F_2nd(:,:,iSym,:) = OF(n).data.F(t2,:,iSym,:);
        OF(n).data.F_1st_Raw(:,:,:,:) = OF(n).data.F(t1,:,:,:);
        OF(n).data.F_2nd_Raw(:,:,:,:) = OF(n).data.F(t2,:,:,:);
        OF(n).data.F_BL(:,:,:,:) = OF(n).data.F(t0,:,:,:); % baseline
        OF(n).data.DFF_1st(:,:,iSym,:) = OF(n).data.DFF(t1,:,iSym,:);
        OF(n).data.DFF_2nd(:,:,iSym,:) = OF(n).data.DFF(t2,:,iSym,:);
        OF(n).data.DFF_1st_Raw(:,:,:,:) = OF(n).data.DFF(t1,:,:,:);
        OF(n).data.DFF_2nd_Raw(:,:,:,:) = OF(n).data.DFF(t2,:,:,:);
        OF(n).data.DFF_BL(:,:,:,:) = OF(n).data.DFF(t0,:,:,:); % baseline
        StimuliFlipped(iSym,1) = Stimuli(iSym,1);
        EpochsFlipped(iSym,1) = 1;

        OF(n).data.F_1st_R(:,:,iSym,:) = OF(n).data.F(t1,:,iSym,:);
        OF(n).data.F_2nd_R(:,:,iSym,:) = OF(n).data.F(t2,:,iSym,:);
        OF(n).data.F_BL_R(:,:,:,:) = OF(n).data.F(t0,:,:,:); % baseline
        OF(n).data.DFF_1st_R(:,:,iSym,:) = OF(n).data.DFF(t1,:,iSym,:);
        OF(n).data.DFF_2nd_R(:,:,iSym,:) = OF(n).data.DFF(t2,:,iSym,:);
        OF(n).data.DFF_BL_R(:,:,:,:) = OF(n).data.DFF(t0,:,:,:); % baseline

        % Pass the data of midline asymmetric OF stimuli
        switch OF(n).somaSide % cell body location of the recorded neuron
            case 'L' % recordings are NOT flipped if soma is on the Left
                OF(n).data.F_1st(:,:,iAsym,:)=...
                    OF(n).data.F(t1,:,iAsym,:);
                OF(n).data.F_2nd(:,:,iAsym,:)=...
                    OF(n).data.F(t2,:,iAsym,:);
                OF(n).data.F_1st(:,:,iProgSide45,:)=...
                    OF(n).data.F(t1,:,iProgSide45,:);
                OF(n).data.F_2nd(:,:,iProgSide45,:)=...
                    OF(n).data.F(t2,:,iProgSide45,:);
                OF(n).data.F_1st(:,:,iProgSide315,:)=...
                    OF(n).data.F(t1,:,iProgSide315,:);
                OF(n).data.F_2nd(:,:,iProgSide315,:)=...
                    OF(n).data.F(t2,:,iProgSide315,:);

                OF(n).data.DFF_1st(:,:,iAsym,:)=...
                    OF(n).data.DFF(t1,:,iAsym,:);
                OF(n).data.DFF_2nd(:,:,iAsym,:)=...
                    OF(n).data.DFF(t2,:,iAsym,:);
                OF(n).data.DFF_1st(:,:,iProgSide45,:)=...
                    OF(n).data.DFF(t1,:,iProgSide45,:);
                OF(n).data.DFF_2nd(:,:,iProgSide45,:)=...
                    OF(n).data.DFF(t2,:,iProgSide45,:);
                OF(n).data.DFF_1st(:,:,iProgSide315,:)=...
                    OF(n).data.DFF(t1,:,iProgSide315,:);
                OF(n).data.DFF_2nd(:,:,iProgSide315,:)=...
                    OF(n).data.DFF(t2,:,iProgSide315,:);

                % note that t2 does not exist for Trans+Rot compound OF
                OF(n).data.F_1st(:,:,iTrans,:)=...
                    OF(n).data.F(t1,:,iTrans,:);
                OF(n).data.DFF_1st(:,:,iTrans,:)=...
                    OF(n).data.DFF(t1,:,iTrans,:);

                % % Flip recordings for the _R variables (see case 'R' for
                % explanation on how to flip)
                OF(n).data.F_1st_R(:,:,iAsym,:)=OF(n).data.F(t2,:,iAsym,:);
                OF(n).data.F_2nd_R(:,:,iAsym,:)=OF(n).data.F(t1,:,iAsym,:);
                OF(n).data.DFF_1st_R(:,:,iAsym,:)=...
                    OF(n).data.DFF(t2,:,iAsym,:);
                OF(n).data.DFF_2nd_R(:,:,iAsym,:)=...
                    OF(n).data.DFF(t1,:,iAsym,:);

                tempRight = OF(n).data.F_1st_R(:,:,RightIndices,:);
                tempLeft = OF(n).data.F_1st_R(:,:,LeftIndices,:);
                OF(n).data.F_1st_R(:,:,RightIndices,:) = tempLeft;
                OF(n).data.F_1st_R(:,:,LeftIndices,:) = tempRight;
                tempRight = OF(n).data.F_2nd_R(:,:,RightIndices,:);
                tempLeft = OF(n).data.F_2nd_R(:,:,LeftIndices,:);
                OF(n).data.F_2nd_R(:,:,RightIndices,:) = tempLeft;
                OF(n).data.F_2nd_R(:,:,LeftIndices,:) = tempRight;

                tempRight = OF(n).data.DFF_1st_R(:,:,RightIndices,:);
                tempLeft = OF(n).data.DFF_1st_R(:,:,LeftIndices,:);
                OF(n).data.DFF_1st_R(:,:,RightIndices,:) = tempLeft;
                OF(n).data.DFF_1st_R(:,:,LeftIndices,:) = tempRight;
                tempRight = OF(n).data.DFF_2nd_R(:,:,RightIndices,:);
                tempLeft = OF(n).data.DFF_2nd_R(:,:,LeftIndices,:);
                OF(n).data.DFF_2nd_R(:,:,RightIndices,:) = tempLeft;
                OF(n).data.DFF_2nd_R(:,:,LeftIndices,:) = tempRight;

                OF(n).data.F_1st_R(:,:,iProgSide45,:)=...
                    OF(n).data.F(t1,:,iProgSide315,:);
                OF(n).data.F_2nd_R(:,:,iProgSide45,:)=...
                    OF(n).data.F(t2,:,iProgSide315,:);
                OF(n).data.F_1st_R(:,:,iProgSide315,:)=...
                    OF(n).data.F(t1,:,iProgSide45,:);
                OF(n).data.F_2nd_R(:,:,iProgSide315,:)=...
                    OF(n).data.F(t2,:,iProgSide45,:);

                OF(n).data.DFF_1st_R(:,:,iProgSide45,:)=...
                    OF(n).data.DFF(t1,:,iProgSide315,:);
                OF(n).data.DFF_2nd_R(:,:,iProgSide45,:)=...
                    OF(n).data.DFF(t2,:,iProgSide315,:);
                OF(n).data.DFF_1st_R(:,:,iProgSide315,:)=...
                    OF(n).data.DFF(t1,:,iProgSide45,:);
                OF(n).data.DFF_2nd_R(:,:,iProgSide315,:)=...
                    OF(n).data.DFF(t2,:,iProgSide45,:);

                OF(n).data.F_1st_R(:,:,iTransCW,:)=...
                    OF(n).data.F(t1,:,iTransCCW,:);
                OF(n).data.F_1st_R(:,:,iTransCCW,:)=...
                    OF(n).data.F(t1,:,iTransCW,:);
                OF(n).data.DFF_1st_R(:,:,iTransCW,:)=...
                    OF(n).data.DFF(t1,:,iTransCCW,:);
                OF(n).data.DFF_1st_R(:,:,iTransCCW,:)=...
                    OF(n).data.DFF(t1,:,iTransCW,:);

            case 'R' % recordings are flipped if soma is on the Right.
                % For midline asymmetric, binocular, rotational OF stimuli
                % - e.g. Yaw OF stimulus moving rightward on 1st epoch and
                % leftward on 2nd epoch - generates front-to-back motion to
                % a Left neuron on the 1st epoch and back-to-front motion
                % on the 2nd epoch. This is the opposite for a Right neuron
                % Therefore, we flip the epochs of Right neurons to match
                % the responses of Left neurons. Use OF.data.DFF_1st_Raw
                % and OF.data.DFF_2nd_Raw to work with unflipped traces.
                OF(n).data.F_1st(:,:,iAsym,:) = OF(n).data.F(t2,:,iAsym,:);
                OF(n).data.F_2nd(:,:,iAsym,:) = OF(n).data.F(t1,:,iAsym,:);

                OF(n).data.DFF_1st(:,:,iAsym,:)=...
                    OF(n).data.DFF(t2,:,iAsym,:);
                OF(n).data.DFF_2nd(:,:,iAsym,:)=...
                    OF(n).data.DFF(t1,:,iAsym,:);
                EpochsFlipped(iAsym,1) = 2;
                StimuliFlipped(iAsym,1) = Stimuli(iAsym,1);

                % If an OF pattern was presented monocularly, (e.g. Yaw OF
                % presented only on Left half of the LED arena) that
                % corresponds to an ipsilateral OF stimulus for a Left
                % neuron and a contralateral OF stimulus for a Right neuron
                % As such, we flip the monocular data for Right neurons,
                % effectively converting them into Left neurons
                tempRight = OF(n).data.F_1st(:,:,RightIndices,:);
                tempLeft = OF(n).data.F_1st(:,:,LeftIndices,:);
                OF(n).data.F_1st(:,:,RightIndices,:) = tempLeft;
                OF(n).data.F_1st(:,:,LeftIndices,:) = tempRight;
                tempRight = OF(n).data.F_2nd(:,:,RightIndices,:);
                tempLeft = OF(n).data.F_2nd(:,:,LeftIndices,:);
                OF(n).data.F_2nd(:,:,RightIndices,:) = tempLeft;
                OF(n).data.F_2nd(:,:,LeftIndices,:) = tempRight;

                tempRight = OF(n).data.DFF_1st(:,:,RightIndices,:);
                tempLeft = OF(n).data.DFF_1st(:,:,LeftIndices,:);
                OF(n).data.DFF_1st(:,:,RightIndices,:) = tempLeft;
                OF(n).data.DFF_1st(:,:,LeftIndices,:) = tempRight;
                tempRight = OF(n).data.DFF_2nd(:,:,RightIndices,:);
                tempLeft = OF(n).data.DFF_2nd(:,:,LeftIndices,:);
                OF(n).data.DFF_2nd(:,:,RightIndices,:) = tempLeft;
                OF(n).data.DFF_2nd(:,:,LeftIndices,:) = tempRight;

                StimuliFlipped(RightIndices,1) = Stimuli(LeftIndices,1);
                StimuliFlipped(LeftIndices,1) = Stimuli(RightIndices,1);

                % For binocular, translational OF patterns that have a
                % sideslip component, the stimuli are midline asymmetric.
                % However, a progressive expanding OF with a rightward
                % sideslip component (ProgSide45) is mirror symmetric to a
                % progressive expanding OF with a leftward sideslip
                % component (ProgSide315). Therefore, we can simply flip
                % the order of these recordings in Right neurons to match
                % the responses of Left neurons
                OF(n).data.F_1st(:,:,iProgSide45,:)=...
                    OF(n).data.F(t1,:,iProgSide315,:);
                OF(n).data.F_2nd(:,:,iProgSide45,:)=...
                    OF(n).data.F(t2,:,iProgSide315,:);
                OF(n).data.F_1st(:,:,iProgSide315,:)=...
                    OF(n).data.F(t1,:,iProgSide45,:);
                OF(n).data.F_2nd(:,:,iProgSide315,:)=...
                    OF(n).data.F(t2,:,iProgSide45,:);

                OF(n).data.DFF_1st(:,:,iProgSide45,:)=...
                    OF(n).data.DFF(t1,:,iProgSide315,:);
                OF(n).data.DFF_2nd(:,:,iProgSide45,:)=...
                    OF(n).data.DFF(t2,:,iProgSide315,:);
                OF(n).data.DFF_1st(:,:,iProgSide315,:)=...
                    OF(n).data.DFF(t1,:,iProgSide45,:);
                OF(n).data.DFF_2nd(:,:,iProgSide315,:)=...
                    OF(n).data.DFF(t2,:,iProgSide45,:);

                OF(n).data.F_1st(:,:,iTransCW,:)=...
                    OF(n).data.F(t1,:,iTransCCW,:);
                OF(n).data.F_1st(:,:,iTransCCW,:)=...
                    OF(n).data.F(t1,:,iTransCW,:);
                OF(n).data.DFF_1st(:,:,iTransCW,:)=...
                    OF(n).data.DFF(t1,:,iTransCCW,:);
                OF(n).data.DFF_1st(:,:,iTransCCW,:)=...
                    OF(n).data.DFF(t1,:,iTransCW,:);


                % Don't Flip it for the _R variables
                OF(n).data.F_1st_R(:,:,iAsym,:)=...
                    OF(n).data.F(t1,:,iAsym,:);
                OF(n).data.F_2nd_R(:,:,iAsym,:)=...
                    OF(n).data.F(t2,:,iAsym,:);
                OF(n).data.F_1st_R(:,:,iProgSide45,:)=...
                    OF(n).data.F(t1,:,iProgSide45,:);
                OF(n).data.F_2nd_R(:,:,iProgSide45,:)=...
                    OF(n).data.F(t2,:,iProgSide45,:);
                OF(n).data.F_1st_R(:,:,iProgSide315,:)=...
                    OF(n).data.F(t1,:,iProgSide315,:);
                OF(n).data.F_2nd_R(:,:,iProgSide315,:)=...
                    OF(n).data.F(t2,:,iProgSide315,:);

                OF(n).data.DFF_1st_R(:,:,iAsym,:)=...
                    OF(n).data.DFF(t1,:,iAsym,:);
                OF(n).data.DFF_2nd_R(:,:,iAsym,:)=...
                    OF(n).data.DFF(t2,:,iAsym,:);
                OF(n).data.DFF_1st_R(:,:,iProgSide45,:)=...
                    OF(n).data.DFF(t1,:,iProgSide45,:);
                OF(n).data.DFF_2nd_R(:,:,iProgSide45,:)=...
                    OF(n).data.DFF(t2,:,iProgSide45,:);
                OF(n).data.DFF_1st_R(:,:,iProgSide315,:)=...
                    OF(n).data.DFF(t1,:,iProgSide315,:);
                OF(n).data.DFF_2nd_R(:,:,iProgSide315,:)=...
                    OF(n).data.DFF(t2,:,iProgSide315,:);

                OF(n).data.F_1st_R(:,:,iTrans,:)=...
                    OF(n).data.F(t1,:,iTrans,:);
                OF(n).data.DFF_1st_R(:,:,iTrans,:)=...
                    OF(n).data.DFF(t1,:,iTrans,:);

                EpochsFlipped(iProgSide45,1) = 1;
                EpochsFlipped(iProgSide315,1) = 1;
                StimuliFlipped(iProgSide45,1) = Stimuli(iProgSide315,1);
                StimuliFlipped(iProgSide315,1) = Stimuli(iProgSide45,1);
                EpochsFlipped(iTransCW,1) = 1;
                EpochsFlipped(iTransCCW,1) = 1;
                StimuliFlipped(iTransCW,1) = Stimuli(iTransCCW,1);
                StimuliFlipped(iTransCCW,1) = Stimuli(iTransCW,1);
        end
    end
end
% remove unnecessary variables
clearvars  i_t0_end i_t0_start i_t1_end i_t1_start i_t2_end i_t2_start ...
    iRollBilat iRollUni iSsBilat iSymBilat iYawBilat iYawUni idx iYF ...
    iSym iSymUni iAsym iLiftBilat iProgSide315 iProgSide45 t t0 t1 t2 n ...
    RightIndices LeftIndices TimeStructure AddTime frameNum blTime ...
    tempLeft tempRight

%% ipsi-, contra-, bi-lateral response matrices for speed tuning

% As of 210701 speed tuning protocols have only one direction of motion (CW
% or CCW) for any given stimulus. Therefore when we flip recordings from
% Right hemisphere neurons onto the Left hemisphere we need to also flip
% the stimulus identity from CW(CCW) to CCW(CW), in addition to flipping
% the side of the monocular stimulation. (clockwise Yaw shown only on the
% right side to a Right neuron is equivalent to counter-clockwise Yaw shown
% only on the left side to a Left neuron). For binocular stimuli, only the
% motion direction (CW<->CCW) is flipped.

% Generate a stimulus name vector which will contain all combinations of
% speed tuning stimuli used after flipping Right hemisphere recordings
StimuliSTPooled = StimuliST;

% Separate visual stimuli indices based on laterality
LeftIndicesST = contains(StimuliST,'_Left_');
RightIndicesST = contains(StimuliST,'_Right_');
BilateralIndicesST = ~LeftIndicesST & ~RightIndicesST;

% Indices of Left CW and CCW stimuli (210701 - This is currently Yaw only)
% needs to be updated if unilateral Roll|Pitch|Translational OF are shown!
LCWIndicesST = contains(StimuliST,'_Left_CW');
LCCWIndicesST = contains(StimuliST,'_Left_CCW');
% Indices of Right CW and CCW stimuli (210701 - This is currently Yaw only)
% needs to be updated if unilateral Roll|Pitch|Translational OF are shown!
RCWIndicesST = contains(StimuliST,'_Right_CW');
RCCWIndicesST = contains(StimuliST,'_Right_CCW');

% indices of asymmetric bilateral stimuli (Yaw)
BYCWIndicesST = contains(StimuliST,'q_CW') & ...
    ~contains(StimuliST,'Flipped') & ...
    ~contains(StimuliST,'Roll') & ~contains(StimuliST,'Prog');
BYCCWIndicesST = contains(StimuliST,'q_CCW') & ...
    ~contains(StimuliST,'Flipped') & ...
    ~contains(StimuliST,'Roll') & ~contains(StimuliST,'Prog');

% yaw flipped indices
YFCWIndicesST = contains(StimuliST,'Flipped')&contains(StimuliST,'q_CW');
YFCCWIndicesST = contains(StimuliST,'Flipped')&contains(StimuliST,'q_CCW');
% bilateral roll indices
BRollCWIndicesST = contains(StimuliST,'Roll') & contains(StimuliST,'q_CW');
BRollCCWIndicesST = contains(StimuliST,'Roll')&contains(StimuliST,'q_CCW');
% bilateral progressive expanding indices
BProgCWIndicesST = contains(StimuliST,'Prog') & contains(StimuliST,'q_CW');

% Stimuli names for each subset (B-binocular ST-speed tuning YF-yawFlipped)
StimuliLeftCWST = StimuliST(LCWIndicesST);
StimuliLeftCCWST = StimuliST(LCCWIndicesST);
StimuliRightCWST = StimuliST(RCWIndicesST);
StimuliRightCCWST = StimuliST(RCCWIndicesST);
StimuliBilatCWST = StimuliST(BYCWIndicesST);
StimuliBilatCCWST = StimuliST(BYCCWIndicesST);
StimuliBRollCWST = StimuliST(BRollCWIndicesST);
StimuliBRollCCWST = StimuliST(BRollCCWIndicesST);
StimuliYFCWST = StimuliST(YFCWIndicesST);
StimuliYFCCWST = StimuliST(YFCCWIndicesST);
StimuliBProgCWST = StimuliST(BProgCWIndicesST);

% Add all the missing stimuli protocols that would result from flipping
% existing speed tuning stimuli presented. This is done even if the
% particular stimulus hasn't been presented to the fly.
if isempty(StimuliLeftCWST) && ~isempty(StimuliRightCCWST)
    StimuliSTPooled = [StimuliSTPooled;...
        strrep(StimuliRightCCWST,'_Right_CCW_','_Left_CW_')];
end
if isempty(StimuliLeftCCWST) && ~isempty(StimuliRightCWST)
    StimuliSTPooled = [StimuliSTPooled;...
        strrep(StimuliRightCWST,'_Right_CW_','_Left_CCW_')];
end
if isempty(StimuliRightCWST) && ~isempty(StimuliLeftCCWST)
    StimuliSTPooled = [StimuliSTPooled;...
        strrep(StimuliLeftCCWST,'_Left_CCW_','_Right_CW_')];
end
if isempty(StimuliRightCCWST) && ~isempty(StimuliLeftCWST)
    StimuliSTPooled = [StimuliSTPooled;...
        strrep(StimuliLeftCWST,'_Left_CW_','_Right_CCW_')];
end
if isempty(StimuliBilatCWST) && ~isempty(StimuliBilatCCWST)
    StimuliSTPooled = [StimuliSTPooled;...
        strrep(StimuliBilatCCWST,'_CCW_','_CW_')];
end
if isempty(StimuliBilatCCWST) && ~isempty(StimuliBilatCWST)
    StimuliSTPooled = [StimuliSTPooled;...
        strrep(StimuliBilatCWST,'_CW_','_CCW_')];
end
if isempty(StimuliBRollCWST) && ~isempty(StimuliBRollCCWST)
    StimuliSTPooled = [StimuliSTPooled;...
        strrep(StimuliBRollCCWST,'_CCW_','_CW_')];
end
if isempty(StimuliBRollCCWST) && ~isempty(StimuliBRollCWST)
    StimuliSTPooled = [StimuliSTPooled;...
        strrep(StimuliBRollCWST,'_CW_','_CCW_')];
end
if isempty(StimuliYFCWST) && ~isempty(StimuliYFCCWST)
    StimuliSTPooled = [StimuliSTPooled;...
        strrep(StimuliYFCCWST,'_CCW_','_CW_')];
end
if isempty(StimuliYFCCWST) && ~isempty(StimuliYFCWST)
    StimuliSTPooled = [StimuliSTPooled;...
        strrep(StimuliYFCWST,'_CW_','_CW_')];
end

% % Get the indices for each subset in the new StimuliSTPooled cell array
% Indices of Left and Right CW and CCW stimuli
LCWIndicesSTP = contains(StimuliSTPooled,'_Left_CW');
LCCWIndicesSTP = contains(StimuliSTPooled,'_Left_CCW');
RCWIndicesSTP = contains(StimuliSTPooled,'_Right_CW');
RCCWIndicesSTP = contains(StimuliSTPooled,'_Right_CCW');
% indices of asymmetric bilateral stimuli (Yaw)
BCWIndicesSTP = contains(StimuliSTPooled,'q_CW') & ...
    ~contains(StimuliSTPooled,'Flipped') & ...
    ~contains(StimuliSTPooled,'Roll') & ~contains(StimuliSTPooled,'Prog');
BCCWIndicesSTP = contains(StimuliSTPooled,'q_CCW') & ...
    ~contains(StimuliSTPooled,'Flipped') & ...
    ~contains(StimuliSTPooled,'Roll') & ~contains(StimuliSTPooled,'Prog');
% yaw flipped indices
YFCWIndicesSTP = contains(StimuliSTPooled,'Flipped') &...
    contains(StimuliSTPooled,'q_CW');
YFCCWIndicesSTP = contains(StimuliSTPooled,'Flipped') &...
    contains(StimuliSTPooled,'q_CCW');
% bilateral roll indices
BRollCWIndicesSTP = contains(StimuliSTPooled,'Roll') &...
    contains(StimuliSTPooled,'q_CW');
BRollCCWIndicesSTP = contains(StimuliSTPooled,'Roll') & ...
    contains(StimuliSTPooled,'q_CCW');
% bilateral progressive expanding indices
BProgCWIndicesSTP = contains(StimuliSTPooled,'Prog') &...
    contains(StimuliSTPooled,'q_CW');

% populate F and DFF matrices with baseline and movement epoch time series
for n = 1:length(OF) % for each neurite
    if ~isempty(OF(n).data.DFFST) % if there is speed tuning imaging data
        % generate the new data matrix that will contain the data in pooled
        % coordinates (where the data is flipped if the cell body is on the
        % right hemisphere)
        OF(n).data.FSTPooled = nan(size(OF(n).data.FST,1),...
            size(OF(n).data.FST,2),length(StimuliSTPooled),...
            size(OF(n).data.FST,4));

        OF(n).data.DFFSTPooled = nan(size(OF(n).data.DFFST,1),...
            size(OF(n).data.DFFST,2),length(StimuliSTPooled),...
            size(OF(n).data.DFFST,4));
        % Data matrix where Right Neuron recordings are NOT flipped
        OF(n).data.FST_Raw = nan(size(OF(n).data.FST,1),...
            size(OF(n).data.FST,2),length(StimuliSTPooled),...
            size(OF(n).data.FST,4));
        OF(n).data.FST_Raw(:,:,contains(StimuliSTPooled,...
            StimuliST),:) = OF(n).data.FST;

        OF(n).data.DFFST_Raw = nan(size(OF(n).data.DFFST,1),...
            size(OF(n).data.DFFST,2),length(StimuliSTPooled),...
            size(OF(n).data.DFFST,4));
        OF(n).data.DFFST_Raw(:,:,contains(StimuliSTPooled,...
            StimuliST),:) = OF(n).data.DFFST;

        switch OF(n).somaSide % cell body location of the recorded neuron
            case 'L' % if the recording was done on a left hand side
                % neuron, no flipping is necessary
                OF(n).data.FSTPooled(:,:,contains(StimuliSTPooled,...
                    StimuliST),:) = OF(n).data.FST;
                OF(n).data.DFFSTPooled(:,:,contains(StimuliSTPooled,...
                    StimuliST),:) = OF(n).data.DFFST;
            case 'R' % if the recording was done on a right hand side
                % neuron, flip all the stimuli accordingly. See the comment in
                % the beginning of this section for detailed info about how the
                % flipping is done
                if sum(RCCWIndicesST)
                    OF(n).data.FSTPooled(:,:,LCWIndicesSTP,:) = ...
                        OF(n).data.FST(:,:,RCCWIndicesST,:);
                    OF(n).data.DFFSTPooled(:,:,LCWIndicesSTP,:) = ...
                        OF(n).data.DFFST(:,:,RCCWIndicesST,:);
                end
                if sum(RCWIndicesST)
                    OF(n).data.FSTPooled(:,:,LCCWIndicesSTP,:) = ...
                        OF(n).data.FST(:,:,RCWIndicesST,:);
                    OF(n).data.DFFSTPooled(:,:,LCCWIndicesSTP,:) = ...
                        OF(n).data.DFFST(:,:,RCWIndicesST,:);
                end
                if sum(LCCWIndicesST)
                    OF(n).data.FSTPooled(:,:,RCWIndicesSTP,:) = ...
                        OF(n).data.FST(:,:,LCCWIndicesST,:);
                    OF(n).data.DFFSTPooled(:,:,RCWIndicesSTP,:) = ...
                        OF(n).data.DFFST(:,:,LCCWIndicesST,:);
                end
                if sum(LCWIndicesST)
                    OF(n).data.FSTPooled(:,:,RCCWIndicesSTP,:) = ...
                        OF(n).data.FST(:,:,LCWIndicesST,:);
                    OF(n).data.DFFSTPooled(:,:,RCCWIndicesSTP,:) = ...
                        OF(n).data.DFFST(:,:,LCWIndicesST,:);
                end
                if sum(BYCCWIndicesST)
                    OF(n).data.FSTPooled(:,:,BCWIndicesSTP,:) = ...
                        OF(n).data.FST(:,:,BYCCWIndicesST,:);
                    OF(n).data.DFFSTPooled(:,:,BCWIndicesSTP,:) = ...
                        OF(n).data.DFFST(:,:,BYCCWIndicesST,:);
                end
                if sum(BYCWIndicesST)
                    OF(n).data.FSTPooled(:,:,BCCWIndicesSTP,:) = ...
                        OF(n).data.FST(:,:,BYCWIndicesST,:);
                    OF(n).data.DFFSTPooled(:,:,BCCWIndicesSTP,:) = ...
                        OF(n).data.DFFST(:,:,BYCWIndicesST,:);
                end
                if sum(BRollCCWIndicesST)
                    OF(n).data.FSTPooled(:,:,BRollCWIndicesSTP,:) = ...
                        OF(n).data.FST(:,:,BRollCCWIndicesST,:);
                    OF(n).data.DFFSTPooled(:,:,BRollCWIndicesSTP,:) = ...
                        OF(n).data.DFFST(:,:,BRollCCWIndicesST,:);
                end
                if sum(BRollCWIndicesST)
                    OF(n).data.FSTPooled(:,:,BRollCCWIndicesSTP,:) = ...
                        OF(n).data.FST(:,:,BRollCWIndicesST,:);
                    OF(n).data.DFFSTPooled(:,:,BRollCCWIndicesSTP,:) = ...
                        OF(n).data.DFFST(:,:,BRollCWIndicesST,:);
                end
                if sum(YFCWIndicesST)
                    OF(n).data.FSTPooled(:,:,YFCWIndicesSTP,:) = ...
                        OF(n).data.FST(:,:,YFCWIndicesST,:);
                    OF(n).data.DFFSTPooled(:,:,YFCWIndicesSTP,:) = ...
                        OF(n).data.DFFST(:,:,YFCWIndicesST,:);
                end
                if sum(YFCCWIndicesST)
                    OF(n).data.FSTPooled(:,:,YFCCWIndicesSTP,:) = ...
                        OF(n).data.FST(:,:,YFCCWIndicesST,:);
                    OF(n).data.DFFSTPooled(:,:,YFCCWIndicesSTP,:) = ...
                        OF(n).data.DFFST(:,:,YFCCWIndicesST,:);
                end
                if sum(BProgCWIndicesST)
                    OF(n).data.FSTPooled(:,:,BProgCWIndicesSTP,:) = ...
                        OF(n).data.FST(:,:,BProgCWIndicesST,:);
                    OF(n).data.DFFSTPooled(:,:,BProgCWIndicesSTP,:) = ...
                        OF(n).data.DFFST(:,:,BProgCWIndicesST,:);
                end
            otherwise
                error(['Soma side not correctly specified for ' ...
                    OF(n).barcode])
        end
    end
end
% remove unnecessary variables
clearvars  LeftIndicesST RightIndicesST BilateralIndicesST LCWIndicesST ...
    LCCWIndicesST RCWIndicesST RCCWIndicesST BCWIndicesST BCCWIndicesST ...
    YFCWIndicesST YFCCWIndicesST BRollCWIndicesST BRollCCWIndicesST ...
    StimuliLeftCWST StimuliLeftCCWST StimuliRightCWST StimuliRightCCWST ...
    StimuliBilatCWST StimuliBilatCCWST StimuliYFCWST StimuliYFCCWST ...
    StimuliBRollCWST StimuliBRollCCWST LCWIndicesSTP LCCWIndicesSTP ...
    RCWIndicesSTP RCCWIndicesSTP BCWIndicesSTP BCCWIndicesSTP n ...
    YFCWIndicesSTP YFCCWIndicesSTP BRollCWIndicesSTP BRollCCWIndicesSTP ...
    BProgCWIndicesSTP BProgCWIndicesST

%% Separate data per neuron, normalize and calculate auc per OF

% In this section all the recordings are separated per cell type.
% Additionally the fluorescence traces are normalized either by maximum
% response or Z-scored and the response amplitude is calculated as auc.

% imaged regions in this dataset (axon or dendrite from each cell type)
regions = unique(strcat({OF.neuron},{OF.region}))';

% Color map for each neuron/region.
cReg = [0.5273 0.5195 0;... % olive (DNp15 > RdlFLPStop +/+)
    0.273 0.195 0;... % dark olive (DNp15 > RdlFLPStop -/-)
    0.5273 0.5195 0;... % olive (DNp15)
    0 0 1;... % blue (H2)
    0.7109 0 0.25;...% red (H2rn)
    0.9 0 0.9;... % magenta (HS > Empty Split-TetxLC)
    0.36 0.05 0.42;... % dark magenta (HS > bIPS-TetxLC)
    1 0 1;... % magenta (HS)
    0 0.4 0;... % green (bIPS>RdlFLPStop +/+)
    0 0.05 0;... % dark green (bIPS>RdlFLPStop -/-)
    0 0.4 0;... % green (bIPS)
    1 0.4 0;... % orange (uLPTCrn)
    ];

% % initialize variables
% OF responses per cell type, where Right recordings are flipped
OFRes_bl = cell(length(regions),length(Stimuli)); % DFF baseline
OFRes_1st = cell(length(regions),length(Stimuli)); % DFF 1st epoch
OFRes_2nd = cell(length(regions),length(Stimuli)); % DFF 2nd epoch
OFRes_bl_Z = cell(length(regions),length(Stimuli)); % Z-scored F baseline
OFRes_1st_Z = cell(length(regions),length(Stimuli)); % Z-scored F 1st epoch
OFRes_2nd_Z = cell(length(regions),length(Stimuli)); % Z-scored F 2nd epoch
OFRes_bl_normM = ... % DFF normalized by max response, baseline
    cell(length(regions),length(Stimuli));
OFRes_1st_normM = ... % DFF normalized by max response, 1st epoch
    cell(length(regions),length(Stimuli));
OFRes_2nd_normM = ... % DFF normalized by max response, 2nd epoch
    cell(length(regions),length(Stimuli));
STResP = cell(length(regions),length(StimuliSTPooled)); % DFF full duration
STResP_Z = ... % Z-scored F full duration
    cell(length(regions),length(StimuliSTPooled));
STResP_normM ... % DFF normalized by max response, full duration
    = cell(length(regions),length(StimuliSTPooled));

% OF responses per cell type, where Right recordings are NOT flipped
OFRes_Raw_bl = cell(length(regions),length(Stimuli)); % DFF baseline
OFRes_Raw_1st = cell(length(regions),length(Stimuli)); % DFF 1st epoch
OFRes_Raw_2nd = cell(length(regions),length(Stimuli)); % DFF 2nd epoch
OFRes_Raw_bl_Z = ...
    cell(length(regions),length(Stimuli)); % Z-scored F baseline
OFRes_Raw_1st_Z = ...
    cell(length(regions),length(Stimuli)); % Z-scored F 1st epoch
OFRes_Raw_2nd_Z = ...
    cell(length(regions),length(Stimuli)); % Z-scored F 2nd epoch
OFRes_Raw_bl_normM = ... % DFF normalized by max response, baseline
    cell(length(regions),length(Stimuli));
OFRes_Raw_1st_normM = ... % DFF normalized by max response, 1st epoch
    cell(length(regions),length(Stimuli));
OFRes_Raw_2nd_normM = ... % DFF normalized by max response, 2nd epoch
    cell(length(regions),length(Stimuli));
STRes = cell(length(regions),length(StimuliST)); % DFF full duration
STRes_Z = ... % Z-scored F full duration
    cell(length(regions),length(StimuliST));
STRes_normM = ... % DFF normalized by max response, full duration
    cell(length(regions),length(StimuliST));
STResP_Raw = ... % DFF full duration
    cell(length(regions),length(StimuliSTPooled));
STResP_Raw_Z = ... % Z-scored F full duration
    cell(length(regions),length(StimuliSTPooled));
STResP_Raw_normM = ... % DFF normalized by max response, full duration
    cell(length(regions),length(StimuliSTPooled));

% OF responses per cell type, where Left recordings are flipped
OFRes_bl_R = cell(length(regions),length(Stimuli)); % DFF baseline
OFRes_1st_R = cell(length(regions),length(Stimuli)); % DFF 1st epoch
OFRes_2nd_R = cell(length(regions),length(Stimuli)); % DFF 2nd epoch
OFRes_bl_Z_R = cell(length(regions),length(Stimuli)); % Z-scored F baseline
OFRes_1st_Z_R=cell(length(regions),length(Stimuli)); % Z-scored F 1st epoch
OFRes_2nd_Z_R=cell(length(regions),length(Stimuli)); % Z-scored F 2nd epoch
OFRes_bl_normM_R = ... % DFF normalized by max response, baseline
    cell(length(regions),length(Stimuli));
OFRes_1st_normM_R = ... % DFF normalized by max response, 1st epcoh
    cell(length(regions),length(Stimuli));
OFRes_2nd_normM_R = ... % DFF normalized by max response, 2nd epcoh
    cell(length(regions),length(Stimuli));
STResP_R=cell(length(regions),length(StimuliSTPooled)); % DFF full duration
STResP_Z_R = ... % Z-scored F full duration
    cell(length(regions),length(StimuliSTPooled));
STResP_normM_R ... % DFF normalized by max response, full duration
    = cell(length(regions),length(StimuliSTPooled));


Auc_1st_Z = cell(length(regions),length(Stimuli));
Auc_2nd_Z = cell(length(regions),length(Stimuli));
Auc_bl_Z = cell(length(regions),length(Stimuli));
Auc_1st_bs_Z = cell(length(regions),length(Stimuli)); % baseline subtracted
Auc_2nd_bs_Z = cell(length(regions),length(Stimuli)); % baseline subtracted
Auc_bl_Z = cell(length(regions),length(Stimuli));
Auc_1st = cell(length(regions),length(Stimuli));
Auc_2nd = cell(length(regions),length(Stimuli));
Auc_1st_bs = cell(length(regions),length(Stimuli)); % baseline subtracted
Auc_2nd_bs = cell(length(regions),length(Stimuli)); % baseline subtracted
Auc_bl = cell(length(regions),length(Stimuli));

Auc_1st_Z_R = cell(length(regions),length(Stimuli));
Auc_2nd_Z_R = cell(length(regions),length(Stimuli));
Auc_bl_Z_R = cell(length(regions),length(Stimuli));
Auc_1st_bs_Z_R=cell(length(regions),length(Stimuli)); % baseline subtracted
Auc_2nd_bs_Z_R=cell(length(regions),length(Stimuli)); % baseline subtracted
Auc_bl_Z_R = cell(length(regions),length(Stimuli));
Auc_1st_R = cell(length(regions),length(Stimuli));
Auc_2nd_R = cell(length(regions),length(Stimuli));
Auc_1st_bs_R = cell(length(regions),length(Stimuli)); % baseline subtracted
Auc_2nd_bs_R = cell(length(regions),length(Stimuli)); % baseline subtracted
Auc_bl_R = cell(length(regions),length(Stimuli));


% Metadata for recorded cells/regions. used for selecting subsets of data
% to plot and analyze
PrepType = cell(length(regions),1);
PrepTypePerROI = cell(length(regions),1);
GCaMPType = cell(length(regions),1);
GCaMPTypePerROI = cell(length(regions),1);
SomaSide = cell(length(regions),1);
SomaSidePerROI = cell(length(regions),1);
Genotype = cell(length(regions),1);
GenotypePerROI = cell(length(regions),1);

numNeurites = zeros(length(regions),1);
numROIs = zeros(length(regions),1);
numNeuritesST = zeros(length(regions),1);
numROIsST = zeros(length(regions),1);

blTime = [1 2]; % Time interval where the baseline mu and sigma are
% calculated for z-scoring
mTime = [2 4]; % Time interval where stimulus is moving

% Define the baseline indices
t=blTime(1);
[~,~,idx]=unique(abs(TimeAxis-t));
i_t0_start = find(idx==1);
% end index for 1st epoch
t=blTime(2);
[~,~,idx]=unique(abs(TimeAxis-t));
i_t0_end = find(idx==1);
% Baseline indices
t0 = i_t0_start:i_t0_end;
% Define the visual motion indices
t=mTime(1);
[~,~,idx]=unique(abs(TimeAxis-t));
i_t1_start = find(idx==1);
% end index for 1st epoch
t=mTime(2);
[~,~,idx]=unique(abs(TimeAxis-t));
i_t1_end = find(idx==1);
% Visual motion indices
tmov = i_t1_start:i_t1_end;

% Populate the matrices with data
for cr = 1:length(regions) % for each cell & region
    OFPerNeuron = OF(contains({OF.barcode},regions{cr}));
    % initialize/reset variables per cell
    meanDFFBilatPD = []; meanDFFIpsiPD = []; meanDFFContraPD = [];
    meanDFFBilatND = []; meanDFFIpsiND = []; meanDFFContraND = [];
    meanF1st = []; meanF2nd = []; meanF1stR = []; meanF2ndR = [];
    meanDFF1st = []; meanDFF2nd = []; meanDFF1stR = []; meanDFF2ndR = [];
    meanDFFST = []; meanDFFSTP = []; meanDFFZST = []; meanDFFZSTP = [];
    STcounterNeurite = 0; STcounterROI = 0; pType = {}; gtype = [];
    sSide = {}; gcType = {}; meanDFFRaw1st = []; meanDFFRaw2nd = [];
    tempFRawZ1st = []; tempFRawZ2nd = []; meanDFFbl = []; meanDFFblR = [];
    meanFRaw1st = []; meanFRaw2nd = []; meanFZ1st = []; meanFZ2nd = [];
    meanFZ1stR = []; meanFZ2ndR = []; meanFRawZ2nd = []; meanFbl = [];
    tempFRawZ = []; tempFRawZ2 = []; meanFRawZ1st = []; meanFblR = [];
    meanFZbl = []; tempFZbl = []; meanFZblR = []; tempFZblR = [];
    maxes = []; maxesST = []; maxesSTP = [];

    if ~isempty(OFPerNeuron) % if there is any data with this neuron class
        % initialize the metadata
        numNeurites(cr) = length(OFPerNeuron);
        PrepType{cr} = cell(length(OFPerNeuron),1);
        GCaMPType{cr} = cell(length(OFPerNeuron),1);
        SomaSide{cr} = cell(length(OFPerNeuron),1);
        Genotype{cr} = cell(length(OFPerNeuron),1);
        for n = 1:numNeurites(cr) % for each neurite
            % calculate population mean and std of the entire imaging
            % session per ROI
            mu = nanmean(OFPerNeuron(n).data.DFF_All_concatenated,2);
            sigma = nanstd(OFPerNeuron(n).data.DFF_All_concatenated,0,2);

            if ~isempty(OFPerNeuron(n).data.DFF) %if there is OF panel data
                % save the calcium indicator used in GCaMPType cell array
                if contains(OFPerNeuron(n).genotype,'sytGCaMP7f') || ...
                        contains(OFPerNeuron(n).genotype,'sytGC7f')
                    GCaMPType{cr}{n,1} = 'sytGC7f';
                else
                    if (contains(OFPerNeuron(n).genotype,'GCaMP7f') || ...
                            contains(OFPerNeuron(n).genotype,'GC7f')) ...
                            && ~contains(OFPerNeuron(n).genotype,'syt')
                        GCaMPType{cr}{n,1} = 'GC7f';
                    else
                        if (contains(OFPerNeuron(n).genotype,'GC6f') || ...
                                contains(OFPerNeuron(n).genotype,...
                                'GCaMP6f'))
                            GCaMPType{cr}{n,1} = 'GC6f';
                        else
                            if (contains(OFPerNeuron(n).genotype,'GC8f')...
                                    || contains(OFPerNeuron(n).genotype,...
                                    'GCaMP8f'))
                                GCaMPType{cr}{n,1} = 'GC8f';
                            end
                        end
                    end
                end
                SomaSide{cr}{n,1} = OFPerNeuron(n).somaSide;
                PrepType{cr}{n,1} = OFPerNeuron(n).prepType;
                Genotype{cr}{n,1} = OFPerNeuron(n).genotype;

                % populate mean responses during the 1st and 2nd motion
                % epochs of optic flow stimuli
                tempF=squeeze((nanmean(OFPerNeuron(n).data.F_1st,2)));
                meanF1st=cat(3,meanF1st,tempF);
                tempF2=squeeze((nanmean(OFPerNeuron(n).data.F_2nd,2)));
                meanF2nd=cat(3,meanF2nd,tempF2);
                tempF=squeeze((nanmean(OFPerNeuron(n).data.F_1st_R,2)));
                meanF1stR=cat(3,meanF1stR,tempF);
                tempF2=squeeze((nanmean(OFPerNeuron(n).data.F_2nd_R,2)));
                meanF2ndR=cat(3,meanF2ndR,tempF2);
                tempF=...
                    squeeze((nanmean(OFPerNeuron(n).data.F_1st_Raw,2)));
                meanFRaw1st=cat(3,meanFRaw1st,tempF);
                tempF2=...
                    squeeze((nanmean(OFPerNeuron(n).data.F_2nd_Raw,2)));
                meanFRaw2nd=cat(3,meanFRaw2nd,tempF2);

                tempDFF=squeeze((nanmean(OFPerNeuron(n).data.DFF_1st,2)));
                meanDFF1st=cat(3,meanDFF1st,tempDFF);
                tempDFF2=squeeze((nanmean(OFPerNeuron(n).data.DFF_2nd,2)));
                meanDFF2nd=cat(3,meanDFF2nd,tempDFF2);
                tempDFF=squeeze((nanmean(OFPerNeuron(n).data.DFF_1st_R,2)));
                meanDFF1stR=cat(3,meanDFF1stR,tempDFF);
                tempDFF2=squeeze((nanmean(OFPerNeuron(n).data.DFF_2nd_R,2)));
                meanDFF2ndR=cat(3,meanDFF2ndR,tempDFF2);
                tempDFF=...
                    squeeze((nanmean(OFPerNeuron(n).data.DFF_1st_Raw,2)));
                meanDFFRaw1st=cat(3,meanDFFRaw1st,tempDFF);
                tempDFF2=...
                    squeeze((nanmean(OFPerNeuron(n).data.DFF_2nd_Raw,2)));
                meanDFFRaw2nd=cat(3,meanDFFRaw2nd,tempDFF2);

                % get the max DFF response per neurite. Used to normalize
                % the data with the maximum calcium response per ROI
                maxes = [maxes squeeze(max(max([tempDFF, tempDFF2])))'];

                % populate mean responses during baseline epoch
                tempF=...
                    squeeze((nanmean(OFPerNeuron(n).data.F_BL,2)));
                meanFbl=cat(3,meanFbl,tempF);
                tempDFF=...
                    squeeze((nanmean(OFPerNeuron(n).data.DFF_BL,2)));
                meanDFFbl=cat(3,meanDFFbl,tempDFF);
                tempF=...
                    squeeze((nanmean(OFPerNeuron(n).data.F_BL_R,2)));
                meanFblR=cat(3,meanFblR,tempF);
                tempDFF=...
                    squeeze((nanmean(OFPerNeuron(n).data.DFF_BL_R,2)));
                meanDFFblR=cat(3,meanDFFblR,tempDFF);

                % calculate z-score per ROI
                tempFZ1st = nan(size(OFPerNeuron(n).data.F_1st));
                tempFZ2nd = nan(size(OFPerNeuron(n).data.F_2nd));
                tempFRawZ1st = nan(size(OFPerNeuron(n).data.F_1st_Raw));
                tempFRawZ2nd = nan(size(OFPerNeuron(n).data.F_2nd_Raw));
                tempFZbl = nan(size(OFPerNeuron(n).data.F_BL));
                tempFZ1stR = nan(size(OFPerNeuron(n).data.F_1st_R));
                tempFZ2ndR = nan(size(OFPerNeuron(n).data.F_2nd_R));
                tempFZblR = nan(size(OFPerNeuron(n).data.F_BL_R));
                for ROI = 1:length(mu) % for each ROI
                    for s = 1:length(Stimuli) % for each stimulus
                        tempFZ1st(:,:,s,ROI) = ...
                            (OFPerNeuron(n).data.F_1st(:,:,s,ROI)-...
                            mean(OFPerNeuron(n).data.F_1st(t0,:,s,ROI)))...
                            ./std(OFPerNeuron(n).data.F_1st(t0,:,s,ROI));
                        tempFZ2nd(:,:,s,ROI) = ...
                            (OFPerNeuron(n).data.F_2nd(:,:,s,ROI)-...
                            mean(OFPerNeuron(n).data.F_2nd(t0,:,s,ROI)))...
                            ./std(OFPerNeuron(n).data.F_2nd(t0,:,s,ROI));
                        tempFZ1stR(:,:,s,ROI) = ...
                            (OFPerNeuron(n).data.F_1st_R(:,:,s,ROI)-...
                            mean(OFPerNeuron(n).data.F_1st_R(t0,:,s,ROI)))...
                            ./std(OFPerNeuron(n).data.F_1st_R(t0,:,s,ROI));
                        tempFZ2ndR(:,:,s,ROI) = ...
                            (OFPerNeuron(n).data.F_2nd_R(:,:,s,ROI)-...
                            mean(OFPerNeuron(n).data.F_2nd_R(t0,:,s,ROI)))...
                            ./std(OFPerNeuron(n).data.F_2nd_R(t0,:,s,ROI));
                        tempFRawZ1st(:,:,s,ROI) = ...
                            (OFPerNeuron(n).data.F_1st_Raw(:,:,s,ROI)-...
                            mean(OFPerNeuron(n).data.F_1st_Raw...
                            (t0,:,s,ROI)))./...
                            std(OFPerNeuron(n).data.F_1st_Raw...
                            (t0,:,s,ROI));
                        tempFRawZ2nd(:,:,s,ROI) = ...
                            (OFPerNeuron(n).data.F_2nd_Raw(:,:,s,ROI)-...
                            mean(OFPerNeuron(n).data.F_2nd_Raw...
                            (t0,:,s,ROI)))./...
                            std(OFPerNeuron(n).data.F_2nd_Raw(t0,:,s,ROI));
                        tempFZbl(:,:,s,ROI) = ...
                            (OFPerNeuron(n).data.F_BL(:,:,s,ROI)-...
                            mean(OFPerNeuron(n).data.F_1st(t0,:,s,ROI)))...
                            ./std(OFPerNeuron(n).data.F_1st(t0,:,s,ROI));
                        tempFZblR(:,:,s,ROI) = ...
                            (OFPerNeuron(n).data.F_BL_R(:,:,s,ROI)-...
                            mean(OFPerNeuron(n).data.F_1st_R(t0,:,s,ROI)))...
                            ./std(OFPerNeuron(n).data.F_1st_R(t0,:,s,ROI));
                    end
                    if isempty(pType)
                        pType = {OFPerNeuron(n).prepType};
                        sSide = {OFPerNeuron(n).somaSide};
                        gcType = {GCaMPType{cr}{n,1}};
                        gtype = {Genotype{cr}{n,1}};
                    else
                        pType = [pType;{OFPerNeuron(n).prepType}];
                        sSide = [sSide;{OFPerNeuron(n).somaSide}];
                        gcType = [gcType;{GCaMPType{cr}{n,1}}];
                        gtype = [gtype;{Genotype{cr}{n,1}}];
                    end
                end
                % save the mean z-scored trace for each visual motion epoch
                tempFZ1st = squeeze(nanmean(tempFZ1st,2));
                tempFZ2nd = squeeze(nanmean(tempFZ2nd,2));
                meanFZ1st=cat(3,meanFZ1st,tempFZ1st);
                meanFZ2nd=cat(3,meanFZ2nd,tempFZ2nd);
                tempFZ1stR = squeeze(nanmean(tempFZ1stR,2));
                tempFZ2ndR = squeeze(nanmean(tempFZ2ndR,2));
                meanFZ1stR=cat(3,meanFZ1stR,tempFZ1stR);
                meanFZ2ndR=cat(3,meanFZ2ndR,tempFZ2ndR);
                tempFRawZ1st = squeeze(nanmean(tempFRawZ1st,2));
                tempFRawZ2nd = squeeze(nanmean(tempFRawZ2nd,2));
                meanFRawZ1st=cat(3,meanFRawZ1st,tempFRawZ1st);
                meanFRawZ2nd=cat(3,meanFRawZ2nd,tempFRawZ2nd);
                tempFZbl = squeeze(nanmean(tempFZbl,2));
                meanFZbl=cat(3,meanFZbl,tempFZbl);
                tempFZblR = squeeze(nanmean(tempFZblR,2));
                meanFZblR=cat(3,meanFZblR,tempFZblR);
            end
            if ~isempty(OFPerNeuron(n).data.DFFST)
                if isempty(PrepType{cr}{n,1})
                    PrepType{cr}{n,1} = OFPerNeuron(n).prepType;
                end
                if isempty(SomaSide{cr}{n,1})
                    SomaSide{cr}{n,1} = OFPerNeuron(n).somaSide;
                end
                if isempty(GCaMPType{cr}{n,1})
                    % save the Ca2+ indicator used in GCaMPType cell array
                    if contains(OFPerNeuron(n).genotype,'sytGC7f') || ...
                            contains(OFPerNeuron(n).genotype,'sytGCaMP7f')
                        GCaMPType{cr}{n,1} = 'sytGC7f';
                    else
                        if(contains(OFPerNeuron(n).genotype,'GCaMP7f')||...
                                contains(...
                                OFPerNeuron(n).genotype,'GC7f')) && ...
                                ~contains(OFPerNeuron(n).genotype,'syt')
                            GCaMPType{cr}{n,1} = 'GC7f';
                        else
                            if (contains(OFPerNeuron(n).genotype,'GC6f')...
                                    ||contains(OFPerNeuron(n).genotype,...
                                    'GCaMP6f'))
                                GCaMPType{cr}{n,1} = 'GC6f';
                            else
                                if (contains(OFPerNeuron(n).genotype,...
                                        'GC8f')|| contains(...
                                        OFPerNeuron(n).genotype,'GCaMP8f'))
                                    GCaMPType{cr}{n,1} = 'GC8f';
                                end
                            end
                        end
                    end
                end
                STcounterNeurite=STcounterNeurite+1;
                % calculate mean traced per speed tuning stimulus
                tempDFFST=squeeze((nanmean(OFPerNeuron(n).data.DFFST,2)));
                maxesST=[maxesST squeeze(max(max(tempDFFST)))'];
                tempDFFSTP=squeeze((nanmean(...
                    OFPerNeuron(n).data.DFFSTPooled,2)));
                maxesSTP=[maxesSTP squeeze(max(max(tempDFFSTP)))'];

                % calculate z-score per ROI
                tempFZST = nan(size(OFPerNeuron(n).data.DFFST));
                for ROI = 1:length(mu)
                    for s = 1:length(StimuliST) % for each stimulus
                        tempFZST(:,:,s,ROI) = ...
                            (OFPerNeuron(n).data.DFFST(:,:,s,ROI)-...
                            mean(OFPerNeuron(n).data.DFFST(t0,:,s,ROI)))...
                            ./std(OFPerNeuron(n).data.DFFST(t0,:,s,ROI));
                    end
                end
                if isempty(pType)
                    for ROI = 1:length(mu)
                        if isempty(pType)
                            pType = {OFPerNeuron(n).prepType};
                            sSide = {OFPerNeuron(n).somaSide};
                            gcType = {GCaMPType{cr}{n,1}};
                            gtype = {Genotype{cr}{n,1}};
                        else
                            pType = [pType;{OFPerNeuron(n).prepType}];
                            sSide = [sSide;{OFPerNeuron(n).somaSide}];
                            gcType = [gcType;{GCaMPType{cr}{n,1}}];
                            gtype = [gtype;{Genotype{cr}{n,1}}];
                        end
                    end
                end
                tempFZST = squeeze(nanmean(tempFZST,2));

                % save the raw and z-scored mean traces for each stimulus
                try
                    meanDFFST=cat(3,meanDFFST,tempDFFST);
                    meanDFFZST=cat(3,meanDFFZST,tempFZST);
                catch
                    try
                        meanDFFST=cat(3,meanDFFST,...
                            tempDFFST(1:size(meanDFFST,1),:,:));
                        meanDFFZST=cat(3,meanDFFZST,...
                            tempFZST(1:size(meanDFFZST,1),:,:));
                    catch
                        if size(tempDFFST,1) < size(meanDFFST,1)
                            tempDFFST(end:size(meanDFFST,1),:) = NaN;
                        end
                        meanDFFST=cat(3,meanDFFST,...
                            tempDFFST(1:size(meanDFFST,1),:,:));
                        if size(tempFZST,1) < size(meanDFFZST,1)
                            tempFZST(end:size(meanDFFZST,1),:) = NaN;
                        end
                        meanDFFZST=cat(3,meanDFFZST,...
                            tempFZST(1:size(meanDFFZST,1),:,:));
                    end
                end

                % calculate z-score per ROI for the pooled Speed Tuning
                % where all possible combinations of CW/CCW stimuli are
                % included as a consequence of flipping right hand side
                % recordings
                tempDFFZSTP = nan(size(OFPerNeuron(n).data.DFFSTPooled));
                for ROI = 1:length(mu)
                    for s = 1:length(StimuliSTPooled) % for each stimulus
                        tempDFFZSTP(:,:,s,ROI) = ...
                            (OFPerNeuron(n).data.DFFSTPooled(:,:,s,ROI)-...
                            mean(OFPerNeuron(n).data.DFFSTPooled(...
                            t0,:,s,ROI)))./std(...
                            OFPerNeuron(n).data.DFFSTPooled(t0,:,s,ROI));
                    end
                end
                tempDFFZSTP = squeeze(nanmean(tempDFFZSTP,2));

                % save the raw and z-scored mean traces for each stimulus
                try
                    meanDFFSTP=cat(3,meanDFFSTP,tempDFFSTP);
                    meanDFFZSTP=cat(3,meanDFFZSTP,tempDFFZSTP);
                catch
                    try
                        meanDFFSTP=cat(3,meanDFFSTP,...
                            tempDFFSTP(1:size(meanDFFSTP,1),:,:));
                    catch
                        if size(tempDFFSTP,1) < size(meanDFFSTP,1)
                            tempDFFSTP(end:size(meanDFFSTP,1),:) = NaN;
                        end
                        meanDFFSTP=cat(3,meanDFFSTP,...
                            tempDFFSTP(1:size(meanDFFSTP,1),:,:));
                    end
                    try
                        meanDFFZSTP=cat(3,meanDFFZSTP,...
                            tempDFFZSTP(1:size(meanDFFZSTP,1),:,:));
                    catch
                        if size(tempDFFZSTP,1) < size(meanDFFZSTP,1)
                            tempDFFZSTP(end:size(meanDFFZSTP,1),:) = NaN;
                        end
                        meanDFFZSTP=cat(3,meanDFFZSTP,...
                            tempDFFZSTP(1:size(meanDFFZSTP,1),:,:));
                    end
                end

            end
            if isempty(PrepType{cr}{n,1})
                PrepType{cr}{n,1} = OFPerNeuron(n).prepType;
            end
            if isempty(SomaSide{cr}{n,1})
                SomaSide{cr}{n,1} = OFPerNeuron(n).somaSide;
            end
            if isempty(GCaMPType{cr}{n,1})
                % save the calcium indicator used in GCaMPType cell array
                if contains(OFPerNeuron(n).genotype,'sytGC7f') || ...
                        contains(OFPerNeuron(n).genotype,'sytGCaMP7f')
                    GCaMPType{cr}{n,1} = 'sytGC7f';
                else
                    if(contains(OFPerNeuron(n).genotype,'GCaMP7f')||...
                            contains(OFPerNeuron(n).genotype,'GC7f'))&&...
                            ~contains(OFPerNeuron(n).genotype,'syt')
                        GCaMPType{cr}{n,1} = 'GC7f';
                    else
                        if (contains(OFPerNeuron(n).genotype,'GC6f')||...
                                contains(OFPerNeuron(n).genotype,...
                                'GCaMP6f'))
                            GCaMPType{cr}{n,1} = 'GC6f';
                        else
                            if (contains(OFPerNeuron(n).genotype,'GC8f')...
                                    || contains(OFPerNeuron(n).genotype,...
                                    'GCaMP8f'))
                                GCaMPType{cr}{n,1} = 'GC8f';
                            end
                        end
                    end
                end
            end
        end
        numNeuritesST(cr) = STcounterNeurite;
        numROIs(cr) = size(meanDFFBilatPD,3);
        PrepTypePerROI{cr} = pType;
        GCaMPTypePerROI{cr} = gcType;
        SomaSidePerROI{cr} = sSide;
        GenotypePerROI{cr} = gtype;

        maxDFF = max(max(max([meanDFFBilatPD,meanDFFIpsiPD,...
            meanDFFContraPD,meanDFFBilatND,meanDFFIpsiND,...
            meanDFFContraND])));

        maxDFF2 = max(max(max([meanDFF1st,meanDFF2nd])));

        for s = 1:length(Stimuli)
            OFRes_1st{cr,s} = squeeze(meanDFF1st(:,s,:));
            OFRes_2nd{cr,s} = squeeze(meanDFF2nd(:,s,:));
            OFRes_bl{cr,s} = squeeze(meanDFFbl(:,s,:));
            OFRes_1st_normM{cr,s} = squeeze(meanDFF1st(:,s,:))./maxes;
            OFRes_2nd_normM{cr,s} = squeeze(meanDFF2nd(:,s,:))./maxes;
            OFRes_1st_Z{cr,s} = squeeze(meanFZ1st(:,s,:));
            OFRes_2nd_Z{cr,s} = squeeze(meanFZ2nd(:,s,:));
            OFRes_1st_R{cr,s} = squeeze(meanDFF1stR(:,s,:));
            OFRes_2nd_R{cr,s} = squeeze(meanDFF2ndR(:,s,:));
            OFRes_bl_R{cr,s} = squeeze(meanDFFblR(:,s,:));
            OFRes_1st_normM_R{cr,s} = squeeze(meanDFF1stR(:,s,:))./maxes;
            OFRes_2nd_normM_R{cr,s} = squeeze(meanDFF2ndR(:,s,:))./maxes;
            OFRes_1st_Z_R{cr,s} = squeeze(meanFZ1stR(:,s,:));
            OFRes_2nd_Z_R{cr,s} = squeeze(meanFZ2ndR(:,s,:));
            OFRes_Raw_1st{cr,s} = squeeze(meanDFFRaw1st(:,s,:));
            OFRes_Raw_2nd{cr,s} = squeeze(meanDFFRaw2nd(:,s,:));
            OFRes_Raw_1st_normM{cr,s} = ...
                squeeze(meanDFFRaw1st(:,s,:))./maxes;
            OFRes_Raw_2nd_normM{cr,s} = ...
                squeeze(meanDFFRaw2nd(:,s,:))./maxes;
            OFRes_Raw_1st_Z{cr,s} = squeeze(meanFRawZ1st(:,s,:));
            OFRes_Raw_2nd_Z{cr,s} = squeeze(meanFRawZ2nd(:,s,:));
            OFRes_bl_Z{cr,s} = squeeze(meanFZbl(:,s,:));
            OFRes_bl_Z_R{cr,s} = squeeze(meanFZblR(:,s,:));

            Auc_1st_Z{cr,s} = trapz(OFRes_1st_Z{cr,s}(tmov,:));
            Auc_2nd_Z{cr,s} = trapz(OFRes_2nd_Z{cr,s}(tmov,:));
            Auc_bl_Z{cr,s} = trapz(OFRes_bl_Z{cr,s});
            Auc_1st{cr,s} = trapz(OFRes_1st{cr,s}(tmov,:));
            Auc_2nd{cr,s} = trapz(OFRes_2nd{cr,s}(tmov,:));
            Auc_bl{cr,s} = trapz(OFRes_bl{cr,s});
            Auc_1st_Z_R{cr,s} = trapz(OFRes_1st_Z_R{cr,s}(tmov,:));
            Auc_2nd_Z_R{cr,s} = trapz(OFRes_2nd_Z_R{cr,s}(tmov,:));
            Auc_bl_Z_R{cr,s} = trapz(OFRes_bl_Z_R{cr,s});
            Auc_1st_R{cr,s} = trapz(OFRes_1st_R{cr,s}(tmov,:));
            Auc_2nd_R{cr,s} = trapz(OFRes_2nd_R{cr,s}(tmov,:));
            Auc_bl_R{cr,s} = trapz(OFRes_bl_R{cr,s});

            Auc_1st_bs_Z{cr,s} = trapz(OFRes_1st_Z{cr,s}(tmov,:)) - ...
                trapz(OFRes_bl_Z{cr,s});
            Auc_2nd_bs_Z{cr,s} = trapz(OFRes_2nd_Z{cr,s}(tmov,:)) - ...
                trapz(OFRes_bl_Z{cr,s});
            Auc_1st_bs{cr,s} = trapz(OFRes_1st{cr,s}(tmov,:)) - ...
                trapz(OFRes_bl{cr,s});
            Auc_2nd_bs{cr,s} = trapz(OFRes_2nd{cr,s}(tmov,:)) - ...
                trapz(OFRes_bl{cr,s});
            Auc_1st_bs_Z_R{cr,s} = trapz(OFRes_1st_Z_R{cr,s}(tmov,:)) - ...
                trapz(OFRes_bl_Z_R{cr,s});
            Auc_2nd_bs_Z_R{cr,s} = trapz(OFRes_2nd_Z_R{cr,s}(tmov,:)) - ...
                trapz(OFRes_bl_Z_R{cr,s});
            Auc_1st_bs_R{cr,s} = trapz(OFRes_1st_R{cr,s}(tmov,:)) - ...
                trapz(OFRes_bl_R{cr,s});
            Auc_2nd_bs_R{cr,s} = trapz(OFRes_2nd_R{cr,s}(tmov,:)) - ...
                trapz(OFRes_bl_R{cr,s});
        end
        for s = 1:size(meanDFFBilatPD,2)
            OFRes_PD{cr,contains(Stimuli,StimuliBinocular{s})}=...
                squeeze(meanDFFBilatPD(:,s,:));
            OFRes_ND{cr,contains(Stimuli,StimuliBinocular{s})}=...
                squeeze(meanDFFBilatND(:,s,:));

            OFRes_PD_normM{cr,contains(Stimuli,StimuliBinocular{s})}=...
                squeeze(meanDFFBilatPD(:,s,:))./maxDFF;
            OFRes_ND_normM{cr,contains(Stimuli,StimuliBinocular{s})}=...
                squeeze(meanDFFBilatND(:,s,:))./maxDFF;

            if s<=size(meanDFFIpsiPD,2)
                OFRes_PD{cr,contains(Stimuli,StimuliLeft{s})}=...
                    squeeze(meanDFFIpsiPD(:,s,:));
                OFRes_PD{cr,contains(Stimuli,StimuliRight{s})}=...
                    squeeze(meanDFFContraPD(:,s,:));
                OFRes_PD_normM{cr,contains(Stimuli,StimuliLeft{s})}=...
                    squeeze(meanDFFIpsiPD(:,s,:))./maxDFF;
                OFRes_PD_normM{cr,contains(Stimuli,StimuliRight{s})}=...
                    squeeze(meanDFFContraPD(:,s,:))./maxDFF;

                OFRes_ND{cr,contains(Stimuli,StimuliLeft{s})}=...
                    squeeze(meanDFFIpsiND(:,s,:));
                OFRes_ND{cr,contains(Stimuli,StimuliRight{s})}=...
                    squeeze(meanDFFContraND(:,s,:));
                OFRes_ND_normM{cr,contains(Stimuli,StimuliLeft{s})}=...
                    squeeze(meanDFFIpsiND(:,s,:))./maxDFF;
                OFRes_ND_normM{cr,contains(Stimuli,StimuliRight{s})}=...
                    squeeze(meanDFFContraND(:,s,:))./maxDFF;
            end
        end
        if ~isempty(meanDFFST)
            numROIsST(cr) = size(meanDFFST,3);
            for s = 1:length(StimuliST)
                STRes{cr,s} = squeeze(meanDFFST(:,s,:));
                STRes_normM{cr,s} = squeeze(meanDFFST(:,s,:))./maxesST;
                STRes_Z{cr,s} = squeeze(meanDFFZST(:,s,:));
            end
            for s = 1:length(StimuliSTPooled)
                STResP{cr,s} = squeeze(meanDFFSTP(:,s,:));
                STResP_normM{cr,s} = squeeze(meanDFFSTP(:,s,:))./maxesSTP;
                STResP_Z{cr,s} = squeeze(meanDFFZSTP(:,s,:));
            end
        end
    end
end


% Initialize normalized z-score matrices
OFRes_1st_Znorm = cell(size(OFRes_1st_Z));
OFRes_2nd_Znorm = cell(size(OFRes_2nd_Z));
OFRes_1st_Znorm_R = cell(size(OFRes_1st_Z_R));
OFRes_2nd_Znorm_R = cell(size(OFRes_2nd_Z_R));
% max z-score values per stimulus per fly
z1Max = cellfun(@nanmax,OFRes_1st_Z,'uni',false);
z2Max = cellfun(@nanmax,OFRes_2nd_Z,'uni',false);
z1MaxR = cellfun(@nanmax,OFRes_1st_Z_R,'uni',false);
z2MaxR = cellfun(@nanmax,OFRes_2nd_Z_R,'uni',false);
for cr = 1:length(regions) % for each neuron / region
    % calculate max z-score value per fly
    zmax1st = nanmax(vertcat(z1Max{cr,:})); % max z-value during 1st epoch
    zmax2nd = nanmax(vertcat(z2Max{cr,:})); % max z-value during 2nd epoch
    zmax = max(zmax1st,zmax2nd);
    zmax1stR = nanmax(vertcat(z1MaxR{cr,:})); % max z-value during 1st epoch
    zmax2ndR = nanmax(vertcat(z2MaxR{cr,:})); % max z-value during 2nd epoch
    zmaxR = max(zmax1stR,zmax2ndR);
    for s = 1:length(Stimuli) % for each stimulus
        % Normalize z-scores based on max value per fly
        OFRes_1st_Znorm{cr,s} = OFRes_1st_Z{cr,s}./zmax;
        OFRes_2nd_Znorm{cr,s} = OFRes_2nd_Z{cr,s}./zmax;
        OFRes_1st_Znorm_R{cr,s} = OFRes_1st_Z_R{cr,s}./zmaxR;
        OFRes_2nd_Znorm_R{cr,s} = OFRes_2nd_Z_R{cr,s}./zmaxR;
    end
end

stimToPlot = find(contains(Stimuli,{'Trans'})&contains(Stimuli,{'5px'}));
stimToPlot(contains(Stimuli(stimToPlot),{'Trans1Rot0'})) = [];
stimOrder = [1,2,3,4,9,10,13,14,17,18,7,8,5,6,11,12,15,16,19,20];% manual reordering based on Trans speed
stimToPlot = stimToPlot(stimOrder);
stimToPlotCW = stimToPlot(find(contains(Stimuli(stimToPlot),{'_CW_'})));
stimToPlotCCW = stimToPlot(find(contains(Stimuli(stimToPlot),{'_CCW_'})));

clearvars cr tempDFF tempDFF2 tempDFFZ tempDFFZ2 tempDFFZST tempDFFZSTP ...
    maxesST maxesSTP maxDFF maxDFF2 meanDFF1st meanDFF2nd ROI ...
    maxes meanDFFBilatND meanDFFBilatPD meanDFFContraND meanDFFIpsiND ...
    meanDFFContraPD meanDFFIpsiPD meanDFFST meanDFFSTP mu s ...
    meanDFFZ1st meanDFFZ2nd meanDFFZST meanDFFZSTP STcounterNeurite n ...
    STcounterROI OFPerNeuron tempDFFST StimuliLeft StimuliRight sigma ...
    tempDFFSTP StimuliBilateral gcType sSide meanDFFRawZ1st ...
    meanDFFRawZ2nd meanDFFRaw1st meanDFFRaw2nd tempDFFRawZ pType gtype ...
    tempDFFRawZ2 z1Max z2Max zmax zmax1st zmax2nd i_t0_start i_t0_end ...
    i_t1_start i_t1_end t z1MaxR z2MaxR zmaxR zmax1stR zmax2ndR idx ...
    meanDFF1stR meanDFF2ndR meanDFFbl meanDFFblR meanF1st meanF1stR ...
    meanF2nd meanF2ndR meanFbl meanFblR meanFRaw1st meanFRaw2nd blTime ...
    meanFRawZ1st meanFRawZ2nd meanFZ1st meanFZ1stR meanFZ2nd meanFZ2ndR ...
    meanFZbl meanFZblR mTime tempF tempF2 tempFRawZ tempFRawZ1st ...
    tempFRawZ2 tempFRawZ2nd tempFZ1st tempFZ1stR tempFZ2nd tempFZ2ndR ...
    tempFZbl tempFZblR tempFZST

%% ----End of pre-processing. The sections below are for plotting data ----

%% Figure 1g
% Here you can specify which calcium imaging traces you'd like to plot for
% visualization purposes. You can select which regions, calcium indicators,
% soma side of the imaged neurons, prep type and OF stimuli (non ST).
% Finally, you can specify if and which normalization should be applied to
% the data

% --------------------------- User Input ----------------------------------
% which OF stimuli to plot. Defined by their index in the Stimuli cell
% array. Plotted in the same order. Leave blank to plot all stimuli.
stimToPlot = [15 14 17 16]; % Fig 1g patterns

% which neurons to plot. Defined by their index in the regions cell array.
% Plotted in the same order. Leave blank to plot all neurons.
neuronsToPlot = [8,3]; % Fig 1g neurons

% Which calcium indicator(s) to plot. Defined by the genotype field in OF
% struct array. Leave blank to plot everything. Note that it is possible to
% pool data from multiple indicators (e.g. gcampsToPlot = {GC6f','GC7f'})
% Accepted arguments: {'GC6f','GC7f','sytGC7f','GC8f'}
% Default: gcampsToPlot = {};
gcampsToPlot = {'GC7f','sytGC7f'};

% Which prep types to plot. Defined by the prepType field in OF struct
% array. Leave blank to plot all experiments. Keep in mind that all of the
% OF sensitivity experiments in walking prep were done when the fly was
% stationary and not moving/grooming/struggling etc.
% Accepted arguments: ['walking','immobilized'] Default: prepToPlot = '';
prepToPlot = 'immobilized';

% Neurons from which hemisphere to plot. Defined by the somaSide field in
% OF struct array. Type 'Both' to plot all neurons. In that case, responses
% from Right hemisphere will be flipped to match the Left hemisphere neural
% responses and plotted together.
% Accepted arguments: ['L','R','Both'] Default: sideToPlot = 'Both';
sideToPlot = 'Both';

% specify normalization type. Type 'DFF' for (F-Fo)/Fo normalization.
% 'Max' will normalize the data per ROI to the maximum DFF response
% recorded in that ROI (see variable maxes).
% 'Zscore' will z-score the data per ROI based on the baseline mean and
% standard deviation for each trial done on that ROI.
% 'Znorm' normalizes Z-score data by the max Z-score recorded from each ROI
% Accepted arguments: ['DFF','Max','Zscore','Znorm']
normType = 'Zscore';

% Specify plot type. Type 'Line' for plotting each fly with a different
% colored line together with the mean trace. 'Shaded' will plot the std as
% a shaded area around the mean trace instead.
% Accepted arguments: ['Line','Shaded']
plotType = 'Shaded';

% Define Y-axis limits for each normalization type.
switch normType
    case 'Max'
        limsY = [-0.5 1.1]; % Y-axis limits
    case 'DFF'
        limsY = [-3 10]; % Y-axis limits
    case 'Zscore'
        switch plotType
            case 'Line'
                limsY = [-6 50]; % Y-axis limits
            case 'Shaded'
                limsY = [-3.5 28]; % Y-axis limits
        end
    case 'Znorm'
        limsY = [-0.4 1.05]; % Y-axis limits
end
%----------------------------end of user input-----------------------------

H2flip = 1;
H2i = find(contains(regions,'H2axon')); % index of H2

% Define the subset of data that will be plotted. This will be used in the
% plotSelectDFF function to plot only the specified subset of data as well
% as to calculate and display number of flies & ROIs being plotted
if isempty(stimToPlot)
    stimToPlot = 1:length(Stimuli);
end
if isempty(neuronsToPlot)
    neuronsToPlot = 1:length(regions);
end

prepIndices = cell(length(regions),1);
gcIndices = cell(length(regions),1);
sideIndices = cell(length(regions),1);
prepIndicesFly = cell(length(regions),1);
gcIndicesFly = cell(length(regions),1);
sideIndicesFly = cell(length(regions),1);

if ~isempty(prepToPlot)
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        if isempty(prepIndices{r})
            prepIndices{r} = ...
                contains(PrepTypePerROI{r},prepToPlot);
            prepIndicesFly{r} = ...
                contains(PrepType{r},prepToPlot);
        else
            prepIndices{r} = (prepIndices{r} | ...
                contains(PrepTypePerROI{r},prepToPlot));
            prepIndicesFly{r} = (prepIndicesFly{r} | ...
                contains(PrepType{r},prepToPlot));
        end
    end
else
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        prepIndices{r} = true(size(PrepTypePerROI{r}));
        prepIndicesFly{r} = true(size(PrepType{r}));
    end
end
if ~isempty(gcampsToPlot)
    for p = 1:length(gcampsToPlot) %for each selected GCaMP type
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            if isempty(gcIndices{r})
                gcIndices{r} = ...
                    ismember(GCaMPTypePerROI{r},gcampsToPlot{p});
                gcIndicesFly{r} = ...
                    ismember(GCaMPType{r},gcampsToPlot{p});
            else
                gcIndices{r} = (gcIndices{r} | ...
                    ismember(GCaMPTypePerROI{r},gcampsToPlot{p}));
                gcIndicesFly{r} = (gcIndicesFly{r} | ...
                    ismember(GCaMPType{r},gcampsToPlot{p}));
            end
        end
    end
else
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        gcIndices{r} = true(size(GCaMPTypePerROI{r}));
        gcIndicesFly{r} = true(size(GCaMPType{r}));
    end
end

switch sideToPlot
    case 'Both'
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            sideIndices{r} = true(size(SomaSidePerROI{r}));
            sideIndicesFly{r} = true(size(SomaSide{r}));
        end
        switch normType
            case 'Zscore'
                dataToPlot1st = OFRes_1st_Z_R;
                dataToPlot2nd = OFRes_2nd_Z_R;
                if H2flip
                    dataToPlot1st(H2i,:) = OFRes_1st_Z(H2i,:);
                    dataToPlot2nd(H2i,:) = OFRes_2nd_Z(H2i,:);
                end
            case 'Max'
                dataToPlot1st = OFRes_1st_normM;
                dataToPlot2nd = OFRes_2nd_normM;
            case 'DFF'
                dataToPlot1st = OFRes_1st;
                dataToPlot2nd = OFRes_2nd;
            case 'Znorm'
                dataToPlot1st = OFRes_1st_Znorm;
                dataToPlot2nd = OFRes_2nd_Znorm;
            otherwise
                error(['Normalization type not specified correctly. ',...
                    'Use normType=[''Max'',''Zscore'',''Znorm'' ',...
                    'or ''Raw'']'])
        end
    case {'L','R'}
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            if isempty(sideIndices{r})
                sideIndices{r} = ...
                    contains(SomaSidePerROI{r},sideToPlot);
                sideIndicesFly{r} = ...
                    contains(SomaSide{r},sideToPlot);
            else
                sideIndices{r} = (sideIndices{r} | ...
                    contains(SomaSidePerROI{r},sideToPlot));
                sideIndicesFly{r} = (sideIndicesFly{r} | ...
                    contains(SomaSide{r},sideToPlot));
            end
        end
        switch normType
            case 'Zscore'
                dataToPlot1st = OFRes_Raw_1st_Z;
                dataToPlot2nd = OFRes_Raw_2nd_Z;
            case 'Max'
                dataToPlot1st = OFRes_Raw_1st_normM;
                dataToPlot2nd = OFRes_Raw_2nd_normM;
            case 'Raw'
                dataToPlot1st = OFRes_Raw_1st;
                dataToPlot2nd = OFRes_Raw_2nd;
            case 'Znorm'
                dataToPlot1st = OFRes_1st_Znorm;
                dataToPlot2nd = OFRes_2nd_Znorm;
            otherwise
                error(['Normalization type not specified correctly. ',...
                    'Use normType=[''Max'',''Zscore'',''Znorm'' ',...
                    'or ''Raw'']'])
        end
    otherwise
        error(['Plotting side not specified correctly. ',...
            'Use sideToPlot = [''L'',''R'' or ''Both'']'])
end

f=plotSelectDFF(dataToPlot1st,dataToPlot2nd,neuronsToPlot,...
    stimToPlot,minFrameRate,cReg,prepIndices,prepIndicesFly,...
    gcIndices,gcIndicesFly,sideIndices,sideIndicesFly,Stimuli,...
    regions,limsY,plotType,normType);
%% Figure 1h

% --------------------------- User Input ----------------------------------
% which neurons to plot. Defined by their index in the regions cell array.
% Plotted in the same order. Leave blank to plot all neurons.

neuronsToPlot = [8,3];

% Which calcium indicator(s) to plot. Defined by the genotype field in OF
% struct array. Leave blank to plot everything. Note that it is possible to
% pool data from multiple indicators (e.g. gcampsToPlot = {GC6f','GC7f'})
% Accepted arguments: {'GC6f','GC7f','sytGC7f'} Default: gcampsToPlot = {};
gcampsToPlot = {'sytGC7f','GC7f'};

% Which prep types to plot. Defined by the prepType field in OF struct
% array. Leave blank to plot all experiments.
% Accepted arguments: ['walking','immobilized'] Default: prepToPlot = '';
prepToPlot = 'immobilized';

% Neurons from which hemisphere to plot. Defined by the somaSide field in
% OF struct array. Type 'Both' to plot all neurons. In that case, responses
% from Right hemisphere will be flipped to match the Left hemisphere neural
% responses and plotted together.
% Accepted arguments: ['L','R','Both'] Default: sideToPlot = 'Both';
sideToPlot = 'Both';

% specify normalization type. Type 'Raw' for no normalization.
% 'Max' will normalize the data per ROI to the maximum response recorded in
% that ROI (see variable maxes).
% 'Zscore' will zscore the data per ROI based on the mean and standard
% deviation of the entire recording session done on that ROI (see variables
% mu and sigma). 'Znorm' normalizes Z-score data by the max Z-score
% observed per fly per neuron
% Accepted arguments: ['Zscore','Raw'] Default: 'Raw';
normType = 'Zscore';

% Specify plot type. Type 'Box' for plotting each fly with a different
% colored dot together with the population mean and sd as boxes.
% Currently only 'Box' is supported.
% Accepted arguments: ['Box']
plotType = 'Box';

% Specify grouping. Type 'neuron' for plotting all stimuli for each neuron
% grouped together. 'stimulus' will group all neurons for each stimulus
% Accepted arguments: ['neuron','stimulus']
groupType = 'neuron';

%----------------------------end of user input-----------------------------

peakbin = 500; % time bin for the averaged peak response (ms)
H2Flip = 0;
iH2 = find(contains(regions,'H2axon')); % index of H2
% index for recordings that will come from the right hemisphere
[~,loc] = ismember(StimuliFlipped,Stimuli);

clearvars l lh h stats
stm = {'Yaw_','Progressive','YawRightFlipped','RollInv','Pitch'};
auci = find(contains(Stimuli,stm) & BinocularIndices);
yi = auci(contains(Stimuli(auci),'Yaw_'));
yfi = auci(contains(Stimuli(auci),'YawRightFlipped'));
pri = auci(contains(Stimuli(auci),'Progressive'));
ssi = find(contains(Stimuli,'Sideslip'));
ps45i = find(contains(Stimuli,'ProgSide45'));

% Define the subset of data that will be plotted. This will be used in the
% plotSelectAUC function to plot only the specified subset of data as well
% as to calculate and display number of flies & ROIs being plotted

if isempty(neuronsToPlot)
    neuronsToPlot = 1:length(regions);
end

prepIndices = cell(length(regions),1);
gcIndices = cell(length(regions),1);
sideIndices = cell(length(regions),1);
prepIndicesFly = cell(length(regions),1);
gcIndicesFly = cell(length(regions),1);
sideIndicesFly = cell(length(regions),1);

if ~isempty(prepToPlot)
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        if isempty(prepIndices{r})
            prepIndices{r} = ...
                contains(PrepTypePerROI{r},prepToPlot);
            prepIndicesFly{r} = ...
                contains(PrepType{r},prepToPlot);
        else
            prepIndices{r} = (prepIndices{r} | ...
                contains(PrepTypePerROI{r},prepToPlot));
            prepIndicesFly{r} = (prepIndicesFly{r} | ...
                contains(PrepType{r},prepToPlot));
        end
    end
else
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        prepIndices{r} = true(size(PrepTypePerROI{r}));
        prepIndicesFly{r} = true(size(PrepType{r}));
    end
end
if ~isempty(gcampsToPlot)
    for p = 1:length(gcampsToPlot) %for each selected GCaMP type
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            if isempty(gcIndices{r})
                gcIndices{r} = ...
                    ismember(GCaMPTypePerROI{r},gcampsToPlot{p});
                gcIndicesFly{r} = ...
                    ismember(GCaMPType{r},gcampsToPlot{p});
            else
                gcIndices{r} = (gcIndices{r} | ...
                    ismember(GCaMPTypePerROI{r},gcampsToPlot{p}));
                gcIndicesFly{r} = (gcIndicesFly{r} | ...
                    ismember(GCaMPType{r},gcampsToPlot{p}));
            end
        end
    end
else
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        gcIndices{r} = true(size(GCaMPTypePerROI{r}));
        gcIndicesFly{r} = true(size(GCaMPType{r}));
    end
end
switch sideToPlot
    case 'Both'
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            sideIndices{r} = true(size(SomaSidePerROI{r}));
            sideIndicesFly{r} = true(size(SomaSide{r}));
        end
        switch normType
            case 'Zscore'
                dataToPlot1st = OFRes_1st_Z;
                dataToPlot2nd = OFRes_2nd_Z;
            case 'Max'
                dataToPlot1st = Auc_1st_normM;
                dataToPlot2nd = Auc_2nd_normM;
            case 'Raw'
                dataToPlot1st = Auc_1st_bs;
                dataToPlot2nd = Auc_2nd_bs;
            otherwise
                error(['Normalization type not specified correctly. ',...
                    'Use normType=[''Max'',''Zscore'', ',...
                    'or ''Raw'']'])
        end
    case {'L','R'}
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            if isempty(sideIndices{r})
                sideIndices{r} = ...
                    contains(SomaSidePerROI{r},sideToPlot);
                sideIndicesFly{r} = ...
                    contains(SomaSide{r},sideToPlot);
            else
                sideIndices{r} = (sideIndices{r} | ...
                    contains(SomaSidePerROI{r},sideToPlot));
                sideIndicesFly{r} = (sideIndicesFly{r} | ...
                    contains(SomaSide{r},sideToPlot));
            end
        end
        switch normType
            case 'Zscore'
                dataToPlot1st = Auc_1st_bs_Z;
                dataToPlot2nd = Auc_2nd_bs_Z;
            case 'Max'
                dataToPlot1st = Auc_1st_normM;
                dataToPlot2nd = Auc_2nd_normM;
            case 'Raw'
                dataToPlot1st = Auc_1st_bs;
                dataToPlot2nd = Auc_2nd_bs;
            otherwise
                error(['Normalization type not specified correctly. ',...
                    'Use normType=[''Max'',''Zscore'', ',...
                    'or ''Raw'']'])
        end
    otherwise
        error(['Plotting side not specified correctly. ',...
            'Use sideToPlot = [''L'',''R'' or ''Both'']'])
end

frameNum = size(dataToPlot1st{1,1},1);
TimeAxisPD = linspace(0,frameNum/(minFrameRate/1000),frameNum);
iPeak = find(TimeAxisPD>4-peakbin/1000,1,'first'):...
    find(TimeAxisPD<4,1,'last');

compNames = {'YawFlipped-Yaw'};
ncomps = length(compNames);
sp=0.5;
x = zeros(ncomps,size(neuronsToPlot,2));
trIndex = cell(ncomps,size(neuronsToPlot,2));
f=figure('WindowState','maximized');
limsY = [-1.2 1.7];
my = zeros(1,ncomps);

for n = 1:ncomps
    for k = neuronsToPlot
        m = find(neuronsToPlot==k);
        x(n,m) = (n-1)*length(neuronsToPlot) + 1 + (m-1)*sp;
        iPlot = prepIndices{k} & gcIndices{k} & sideIndices{k};
        switch compNames{n}
            case 'Prog-Yaw'
                switch regions{k}
                    case {'H2axon'}
                        PKprog = mean(dataToPlot2nd{k,pri}(iPeak,iPlot))-mean(dataToPlot1st{k,pri}(iPeak,iPlot));
                        tPKprog = abs(mean(dataToPlot2nd{k,pri}(iPeak,iPlot))-mean(dataToPlot1st{k,pri}(iPeak,iPlot)));
                        PKyaw = mean(dataToPlot1st{k,yi}(iPeak,iPlot))-mean(dataToPlot2nd{k,yi}(iPeak,iPlot));
                        tPKyaw = abs(mean(dataToPlot1st{k,yi}(iPeak,iPlot))-mean(dataToPlot2nd{k,yi}(iPeak,iPlot)));
                        trIndex{n,m}=(PKprog-PKyaw)./(tPKprog+tPKyaw);
                    otherwise
                        PKprog = mean(dataToPlot1st{k,pri}(iPeak,iPlot))-mean(dataToPlot2nd{k,pri}(iPeak,iPlot));
                        tPKprog = abs(mean(dataToPlot1st{k,pri}(iPeak,iPlot))-mean(dataToPlot2nd{k,pri}(iPeak,iPlot)));
                        PKyaw = mean(dataToPlot2nd{k,yi}(iPeak,iPlot))-mean(dataToPlot1st{k,yi}(iPeak,iPlot));
                        tPKyaw = abs(mean(dataToPlot2nd{k,yi}(iPeak,iPlot))-mean(dataToPlot1st{k,yi}(iPeak,iPlot)));
                        trIndex{n,m}=(PKprog-PKyaw)./(tPKprog+tPKyaw);
                end
            case 'YawFlipped-Yaw'
                switch regions{k}
                    case {'H2axon'}
                        PKyf = mean(dataToPlot1st{k,yfi}(iPeak,iPlot))-mean(dataToPlot2nd{k,yfi}(iPeak,iPlot));
                        tPKyf = abs(mean(dataToPlot2nd{k,yfi}(iPeak,iPlot))-mean(dataToPlot1st{k,yfi}(iPeak,iPlot)));
                        PKyaw = mean(dataToPlot1st{k,yi}(iPeak,iPlot))-mean(dataToPlot2nd{k,yi}(iPeak,iPlot));
                        tPKyaw = abs(mean(dataToPlot1st{k,yi}(iPeak,iPlot))-mean(dataToPlot2nd{k,yi}(iPeak,iPlot)));
                        trIndex{n,m}=(PKyf-PKyaw)./(tPKyf+tPKyaw);
                    case {'H2rnaxon'}
                        PKyf = mean(dataToPlot1st{k,yfi}(iPeak,iPlot))-mean(dataToPlot2nd{k,yfi}(iPeak,iPlot));
                        tPKyf = abs(mean(dataToPlot2nd{k,yfi}(iPeak,iPlot))-mean(dataToPlot1st{k,yfi}(iPeak,iPlot)));
                        PKyaw = mean(dataToPlot2nd{k,yi}(iPeak,iPlot))-mean(dataToPlot1st{k,yi}(iPeak,iPlot));
                        tPKyaw = abs(mean(dataToPlot1st{k,yi}(iPeak,iPlot))-mean(dataToPlot2nd{k,yi}(iPeak,iPlot)));
                        trIndex{n,m}=(PKyf-PKyaw)./(tPKyf+tPKyaw);
                    otherwise
                        PKyf = mean(dataToPlot2nd{k,yfi}(iPeak,iPlot))-mean(dataToPlot1st{k,yfi}(iPeak,iPlot));
                        tPKyf = abs(mean(dataToPlot1st{k,yfi}(iPeak,iPlot))-mean(dataToPlot2nd{k,yfi}(iPeak,iPlot)));
                        PKyaw = mean(dataToPlot2nd{k,yi}(iPeak,iPlot))-mean(dataToPlot1st{k,yi}(iPeak,iPlot));
                        tPKyaw = abs(mean(dataToPlot2nd{k,yi}(iPeak,iPlot))-mean(dataToPlot1st{k,yi}(iPeak,iPlot)));
                        trIndex{n,m}=(PKyf-PKyaw)./(tPKyf+tPKyaw);
                end

            case 'Prog-Sideslip'
                switch regions{k}
                    case {'H2axon'}
                        PKss = mean(dataToPlot2nd{k,ssi}(iPeak,iPlot))-mean(dataToPlot1st{k,ssi}(iPeak,iPlot));
                        tPKss = abs(mean(dataToPlot2nd{k,ssi}(iPeak,iPlot))-mean(dataToPlot1st{k,ssi}(iPeak,iPlot)));
                        PKprog = mean(dataToPlot2nd{k,pri}(iPeak,iPlot))-mean(dataToPlot1st{k,pri}(iPeak,iPlot));
                        tPKprog = abs(mean(dataToPlot2nd{k,pri}(iPeak,iPlot))-mean(dataToPlot1st{k,pri}(iPeak,iPlot)));
                        trIndex{n,m}=(PKprog-PKss)./(tPKprog+tPKss);
                    otherwise
                        PKss = mean(dataToPlot1st{k,ssi}(iPeak,iPlot))-mean(dataToPlot2nd{k,ssi}(iPeak,iPlot));
                        tPKss = abs(mean(dataToPlot1st{k,ssi}(iPeak,iPlot))-mean(dataToPlot2nd{k,ssi}(iPeak,iPlot)));
                        PKprog = mean(dataToPlot1st{k,pri}(iPeak,iPlot))-mean(dataToPlot2nd{k,pri}(iPeak,iPlot));
                        tPKprog = abs(mean(dataToPlot1st{k,pri}(iPeak,iPlot))-mean(dataToPlot2nd{k,pri}(iPeak,iPlot)));
                        trIndex{n,m}=(PKprog-PKss)./(tPKprog+tPKss);
                end

            case 'ProgSide45-Sideslip'
                switch regions{k}
                    case {'H2axon'}
                        PKprog45 = mean(dataToPlot2nd{k,ps45i}(iPeak,iPlot))-mean(dataToPlot1st{k,ps45i}(iPeak,iPlot));
                        tPKprog45 = abs(mean(dataToPlot2nd{k,ps45i}(iPeak,iPlot))-mean(dataToPlot1st{k,ps45i}(iPeak,iPlot)));
                        PKss = mean(dataToPlot2nd{k,ssi}(iPeak,iPlot))-mean(dataToPlot1st{k,ssi}(iPeak,iPlot));
                        tPKss = abs(mean(dataToPlot2nd{k,ssi}(iPeak,iPlot))-mean(dataToPlot1st{k,ssi}(iPeak,iPlot)));
                        trIndex{n,m}=(PKprog45-PKss)./(tPKprog45+tPKss);
                    case {'H2rnaxon'}
                        PKprog45 = mean(dataToPlot1st{k,ps45i}(iPeak,iPlot))-mean(dataToPlot2nd{k,ps45i}(iPeak,iPlot));
                        tPKprog45 = abs(mean(dataToPlot1st{k,ps45i}(iPeak,iPlot))-mean(dataToPlot2nd{k,ps45i}(iPeak,iPlot)));
                        PKss = mean(dataToPlot1st{k,ssi}(iPeak,iPlot))-mean(dataToPlot2nd{k,ssi}(iPeak,iPlot));
                        tPKss = abs(mean(dataToPlot1st{k,ssi}(iPeak,iPlot))-mean(dataToPlot2nd{k,ssi}(iPeak,iPlot)));
                        trIndex{n,m}=(PKprog45-PKss)./(tPKprog45+tPKss);
                    otherwise
                        PKprog45 = mean(dataToPlot1st{k,ps45i}(iPeak,iPlot))-mean(dataToPlot2nd{k,ps45i}(iPeak,iPlot));
                        tPKprog45 = abs(mean(dataToPlot1st{k,ps45i}(iPeak,iPlot))-mean(dataToPlot2nd{k,ps45i}(iPeak,iPlot)));
                        PKss = mean(dataToPlot1st{k,ssi}(iPeak,iPlot))-mean(dataToPlot2nd{k,ssi}(iPeak,iPlot));
                        tPKss = abs(mean(dataToPlot1st{k,ssi}(iPeak,iPlot))-mean(dataToPlot2nd{k,ssi}(iPeak,iPlot)));
                        trIndex{n,m}=(PKprog45-PKss)./(tPKprog45+tPKss);
                end
        end
        h(n,m) = notBoxPlot(trIndex{n,m},x(n,m));hold on
        set(h(n,m).data,'MarkerFaceColor',cReg(k,:))
        % save the maximum y values for significance plot alignment
        maxbl = nanmax(trIndex{n,m});
        if isempty(maxbl)
            maxbl = 0;
        end
        my(n) = max(my(n),maxbl);
        trIndex{n,m}(isnan(trIndex{n,m})) = [];
    end
    % Perform Mann-Whitney-Wilcoxon non parametric test and plot
    % significance indicators if p<0.05
    nck = nchoosek(1:m,2);
    %This is all done just to help with the plots
    nck = [nck; nck(end,:)];
    nph = size(nck,1) - 1; % number of comparisons made (used for adjusting the
    % minimum P value required for significance - a.k.a. Bonferroni correction)
    sameComparison = 1;
    for ii = 1 : size(nck,1) - 1
        if nck(ii+1) == nck(ii)
            sameComparison = sameComparison + 1;
        else
            sameComparison = 1;
        end
    end
    for ii = 1 : size(nck,1) - 1
        stats{ii} = mwwtest(trIndex{n,nck(ii,1)},trIndex{n,nck(ii,2)},0);
        % disp([compNames{n},'-',regions{neuronsToPlot(nck(ii,1))},'-',...
        %     regions{neuronsToPlot(nck(ii,2))},'-',num2str(stats{ii}.p(2))])
        if stats{ii}.p < 0.05/nph
            % uncomment below if you want to display the p-values per pair
            disp([compNames{n},'-',regions{neuronsToPlot(nck(ii,1))},'-',...
                regions{neuronsToPlot(nck(ii,2))},'-',num2str(stats{ii}.p(2))])
            sigline([x(n,nck(ii,1)),x(n,nck(ii,2))],[],my(n)+my(n)*0.05*ii);
        end
    end
end
limsX = xlim;
line(limsX,[0,0],'LineStyle','--','Color','k','LineWidth',1)
set(gca,'XTick',x(:,end)','XTickLabel',compNames,'FontSize',18,...
    'YLim',limsY)
set(gca,'FontSize',14,'YLim',limsY)
ylabel('Discrimination Index (peak)')

for m=1:length(neuronsToPlot)
    hold on
    l(m)=plot(nan,nan,'Color',cReg(neuronsToPlot(m),:),'LineWidth',4);
end
set(gca,'XLim',limsX)
lh=legend(l,regions(neuronsToPlot),'Box','off','FontSize',24,'Location',...
    'Best');
clearvars l lh h

%% Figure 3c

% --------------------------- User Input ----------------------------------
% which OF stimuli to plot. Defined by their index in the Stimuli cell
% array. Plotted in the same order. Leave blank to plot all stimuli.

stimToPlot = [15 17 16 14 7]; % Fig 3c patterns

% which neurons to plot. Defined by their index in the regions cell array.
% Plotted in the same order. Leave blank to plot all neurons.

neuronsToPlot = [8,4,11,5,12,3]; % Fig 3c neurons

% Which calcium indicator(s) to plot. Defined by the genotype field in OF
% struct array. Leave blank to plot everything. Note that it is possible to
% pool data from multiple indicators (e.g. gcampsToPlot = {GC6f','GC7f'})
% Accepted arguments: {'GC6f','GC7f','sytGC7f','GC8f'}
% Default: gcampsToPlot = {};
gcampsToPlot = {'GC7f','sytGC7f'};


% Which prep types to plot. Defined by the prepType field in OF struct
% array. Leave blank to plot all experiments. Keep in mind that all of the
% OF sensitivity experiments in walking prep were done when the fly was
% stationary and not moving/grooming/struggling etc.
% Accepted arguments: ['walking','immobilized'] Default: prepToPlot = '';
prepToPlot = 'immobilized';

% Neurons from which hemisphere to plot. Defined by the somaSide field in
% OF struct array. Type 'Both' to plot all neurons. In that case, responses
% from Right hemisphere will be flipped to match the Left hemisphere neural
% responses and plotted together.
% Accepted arguments: ['L','R','Both'] Default: sideToPlot = 'Both';
sideToPlot = 'Both';

% specify normalization type. Type 'DFF' for (F-Fo)/Fo normalization.
% 'Max' will normalize the data per ROI to the maximum DFF response
% recorded in that ROI (see variable maxes).
% 'Zscore' will z-score the data per ROI based on the baseline mean and
% standard deviation for each trial done on that ROI.
% 'Znorm' normalizes Z-score data by the max Z-score recorded from each ROI
% Accepted arguments: ['DFF','Max','Zscore','Znorm']
normType = 'Zscore';

% Specify plot type. Type 'Line' for plotting each fly with a different
% colored line together with the mean trace. 'Shaded' will plot the std as
% a shaded area around the mean trace instead.
% Accepted arguments: ['Line','Shaded']
plotType = 'Shaded';

% Define Y-axis limits for each normalization type.
switch normType
    case 'Max'
        limsY = [-0.5 1.1]; % Y-axis limits
    case 'DFF'
        limsY = [-2 10]; % Y-axis limits
    case 'Zscore'
        switch plotType
            case 'Line'
                limsY = [-6 65]; % Y-axis limits
            case 'Shaded'
                limsY = [-15 32]; % Y-axis limits
        end
    case 'Znorm'
        limsY = [-0.4 1.05]; % Y-axis limits
end
%----------------------------end of user input-----------------------------

H2flip = 1;
H2i = find(contains(regions,'H2axon')); % index of H2

% Define the subset of data that will be plotted. This will be used in the
% plotSelectDFF function to plot only the specified subset of data as well
% as to calculate and display number of flies & ROIs being plotted
if isempty(stimToPlot)
    stimToPlot = 1:length(Stimuli);
end
if isempty(neuronsToPlot)
    neuronsToPlot = 1:length(regions);
end

prepIndices = cell(length(regions),1);
gcIndices = cell(length(regions),1);
sideIndices = cell(length(regions),1);
prepIndicesFly = cell(length(regions),1);
gcIndicesFly = cell(length(regions),1);
sideIndicesFly = cell(length(regions),1);

if ~isempty(prepToPlot)
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        if isempty(prepIndices{r})
            prepIndices{r} = ...
                contains(PrepTypePerROI{r},prepToPlot);
            prepIndicesFly{r} = ...
                contains(PrepType{r},prepToPlot);
        else
            prepIndices{r} = (prepIndices{r} | ...
                contains(PrepTypePerROI{r},prepToPlot));
            prepIndicesFly{r} = (prepIndicesFly{r} | ...
                contains(PrepType{r},prepToPlot));
        end
    end
else
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        prepIndices{r} = true(size(PrepTypePerROI{r}));
        prepIndicesFly{r} = true(size(PrepType{r}));
    end
end
if ~isempty(gcampsToPlot)
    for p = 1:length(gcampsToPlot) %for each selected GCaMP type
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            if isempty(gcIndices{r})
                gcIndices{r} = ...
                    ismember(GCaMPTypePerROI{r},gcampsToPlot{p});
                gcIndicesFly{r} = ...
                    ismember(GCaMPType{r},gcampsToPlot{p});
            else
                gcIndices{r} = (gcIndices{r} | ...
                    ismember(GCaMPTypePerROI{r},gcampsToPlot{p}));
                gcIndicesFly{r} = (gcIndicesFly{r} | ...
                    ismember(GCaMPType{r},gcampsToPlot{p}));
            end
        end
    end
else
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        gcIndices{r} = true(size(GCaMPTypePerROI{r}));
        gcIndicesFly{r} = true(size(GCaMPType{r}));
    end
end

switch sideToPlot
    case 'Both'
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            sideIndices{r} = true(size(SomaSidePerROI{r}));
            sideIndicesFly{r} = true(size(SomaSide{r}));
        end
        switch normType
            case 'Zscore'
                dataToPlot1st = OFRes_1st_Z_R;
                dataToPlot2nd = OFRes_2nd_Z_R;
                if H2flip
                    dataToPlot1st(H2i,:) = OFRes_1st_Z(H2i,:);
                    dataToPlot2nd(H2i,:) = OFRes_2nd_Z(H2i,:);
                end
            case 'Max'
                dataToPlot1st = OFRes_1st_normM;
                dataToPlot2nd = OFRes_2nd_normM;
            case 'DFF'
                dataToPlot1st = OFRes_1st_R;
                dataToPlot2nd = OFRes_2nd_R;
                if H2flip
                    dataToPlot1st(H2i,:) = OFRes_1st(H2i,:);
                    dataToPlot2nd(H2i,:) = OFRes_2nd(H2i,:);
                end
            case 'Znorm'
                dataToPlot1st = OFRes_1st_Znorm;
                dataToPlot2nd = OFRes_2nd_Znorm;
            otherwise
                error(['Normalization type not specified correctly. ',...
                    'Use normType=[''Max'',''Zscore'',''Znorm'' ',...
                    'or ''Raw'']'])
        end
    case {'L','R'}
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            if isempty(sideIndices{r})
                sideIndices{r} = ...
                    contains(SomaSidePerROI{r},sideToPlot);
                sideIndicesFly{r} = ...
                    contains(SomaSide{r},sideToPlot);
            else
                sideIndices{r} = (sideIndices{r} | ...
                    contains(SomaSidePerROI{r},sideToPlot));
                sideIndicesFly{r} = (sideIndicesFly{r} | ...
                    contains(SomaSide{r},sideToPlot));
            end
        end
        switch normType
            case 'Zscore'
                dataToPlot1st = OFRes_Raw_1st_Z;
                dataToPlot2nd = OFRes_Raw_2nd_Z;
            case 'Max'
                dataToPlot1st = OFRes_Raw_1st_normM;
                dataToPlot2nd = OFRes_Raw_2nd_normM;
            case 'Raw'
                dataToPlot1st = OFRes_Raw_1st;
                dataToPlot2nd = OFRes_Raw_2nd;
            case 'Znorm'
                dataToPlot1st = OFRes_1st_Znorm;
                dataToPlot2nd = OFRes_2nd_Znorm;
            otherwise
                error(['Normalization type not specified correctly. ',...
                    'Use normType=[''Max'',''Zscore'',''Znorm'' ',...
                    'or ''Raw'']'])
        end
    otherwise
        error(['Plotting side not specified correctly. ',...
            'Use sideToPlot = [''L'',''R'' or ''Both'']'])
end

f=plotSelectDS(dataToPlot1st,dataToPlot2nd,neuronsToPlot,...
    stimToPlot,minFrameRate,cReg,prepIndices,prepIndicesFly,...
    gcIndices,gcIndicesFly,sideIndices,sideIndicesFly,Stimuli,...
    regions,limsY,plotType,normType);

%% Figure 3d

% --------------------------- User Input ----------------------------------
% which OF stimuli to plot. Defined by their index in the Stimuli cell
% array. Plotted in the same order. Leave blank to plot all stimuli.
stimToPlot = [15 17 16 14 7];

% which neurons to plot. Defined by their index in the regions cell array.
% Plotted in the same order. Leave blank to plot all neurons.

neuronsToPlot = [8,4,11,5,12,3];

% Which calcium indicator(s) to plot. Defined by the genotype field in OF
% struct array. Leave blank to plot everything. Note that it is possible to
% pool data from multiple indicators (e.g. gcampsToPlot = {GC6f','GC7f'})
% Accepted arguments: {'GC6f','GC7f','sytGC7f'} Default: gcampsToPlot = {};
gcampsToPlot = {'sytGC7f','GC7f'};

% Which prep types to plot. Defined by the prepType field in OF struct
% array. Leave blank to plot all experiments.
% Accepted arguments: ['walking','immobilized'] Default: prepToPlot = '';
prepToPlot = 'immobilized';

% Neurons from which hemisphere to plot. Defined by the somaSide field in
% OF struct array. Type 'Both' to plot all neurons. In that case, responses
% from Right hemisphere will be flipped to match the Left hemisphere neural
% responses and plotted together.
% Accepted arguments: ['L','R','Both'] Default: sideToPlot = 'Both';
sideToPlot = 'Both';

% specify normalization type. Type 'Raw' for no normalization.
% 'Max' will normalize the data per ROI to the maximum response recorded in
% that ROI (see variable maxes).
% 'Zscore' will zscore the data per ROI based on the mean and standard
% deviation of the entire recording session done on that ROI (see variables
% mu and sigma). 'Znorm' normalizes Z-score data by the max Z-score
% observed per fly per neuron
% Accepted arguments: ['Zscore','Raw'] Default: 'Raw';
normType = 'Zscore';

%----------------------------end of user input-----------------------------

peakbin = 500; % time bin for the averaged peak response (ms)
iH2 = find(contains(regions,'H2axon')); % index of H2
% index for recordings that will come from the right hemisphere
[~,loc] = ismember(StimuliFlipped,Stimuli);
pairedLines = 1;

clearvars l lh h stats
stm = {'Yaw_','Progressive','YawRightFlipped','RollInv','Pitch'};
auci = find(contains(Stimuli,stm) & BinocularIndices);
yi = auci(contains(Stimuli(auci),'Yaw_')); % index of the stimulus that
% will be used for normalization (currently is binocular yaw)

if isempty(stimToPlot)
    stimToPlot = 1:length(Stimuli);
end
if isempty(neuronsToPlot)
    neuronsToPlot = 1:length(regions);
end

prepIndices = cell(length(regions),1);
gcIndices = cell(length(regions),1);
sideIndices = cell(length(regions),1);
prepIndicesFly = cell(length(regions),1);
gcIndicesFly = cell(length(regions),1);
sideIndicesFly = cell(length(regions),1);

if ~isempty(prepToPlot)
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        if isempty(prepIndices{r})
            prepIndices{r} = ...
                contains(PrepTypePerROI{r},prepToPlot);
            prepIndicesFly{r} = ...
                contains(PrepType{r},prepToPlot);
        else
            prepIndices{r} = (prepIndices{r} | ...
                contains(PrepTypePerROI{r},prepToPlot));
            prepIndicesFly{r} = (prepIndicesFly{r} | ...
                contains(PrepType{r},prepToPlot));
        end
    end
else
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        prepIndices{r} = true(size(PrepTypePerROI{r}));
        prepIndicesFly{r} = true(size(PrepType{r}));
    end
end
if ~isempty(gcampsToPlot)
    for p = 1:length(gcampsToPlot) %for each selected GCaMP type
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            if isempty(gcIndices{r})
                gcIndices{r} = ...
                    ismember(GCaMPTypePerROI{r},gcampsToPlot{p});
                gcIndicesFly{r} = ...
                    ismember(GCaMPType{r},gcampsToPlot{p});
            else
                gcIndices{r} = (gcIndices{r} | ...
                    ismember(GCaMPTypePerROI{r},gcampsToPlot{p}));
                gcIndicesFly{r} = (gcIndicesFly{r} | ...
                    ismember(GCaMPType{r},gcampsToPlot{p}));
            end
        end
    end
else
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        gcIndices{r} = true(size(GCaMPTypePerROI{r}));
        gcIndicesFly{r} = true(size(GCaMPType{r}));
    end
end
switch sideToPlot
    case 'Both'
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            sideIndices{r} = true(size(SomaSidePerROI{r}));
            sideIndicesFly{r} = true(size(SomaSide{r}));
        end
        switch normType
            case 'Zscore'
                dataToPlot1st = OFRes_1st_Z;
                dataToPlot2nd = OFRes_2nd_Z;
            case 'Max'
                dataToPlot1st = Auc_1st_normM;
                dataToPlot2nd = Auc_2nd_normM;
            case 'Raw'
                dataToPlot1st = Auc_1st_bs;
                dataToPlot2nd = Auc_2nd_bs;
            otherwise
                error(['Normalization type not specified correctly. ',...
                    'Use normType=[''Max'',''Zscore'', ',...
                    'or ''Raw'']'])
        end
    case {'L','R'}
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            if isempty(sideIndices{r})
                sideIndices{r} = ...
                    contains(SomaSidePerROI{r},sideToPlot);
                sideIndicesFly{r} = ...
                    contains(SomaSide{r},sideToPlot);
            else
                sideIndices{r} = (sideIndices{r} | ...
                    contains(SomaSidePerROI{r},sideToPlot));
                sideIndicesFly{r} = (sideIndicesFly{r} | ...
                    contains(SomaSide{r},sideToPlot));
            end
        end
        switch normType
            case 'Zscore'
                dataToPlot1st = Auc_1st_bs_Z;
                dataToPlot2nd = Auc_2nd_bs_Z;
            case 'Max'
                dataToPlot1st = Auc_1st_normM;
                dataToPlot2nd = Auc_2nd_normM;
            case 'Raw'
                dataToPlot1st = Auc_1st_bs;
                dataToPlot2nd = Auc_2nd_bs;
            otherwise
                error(['Normalization type not specified correctly. ',...
                    'Use normType=[''Max'',''Zscore'', ',...
                    'or ''Raw'']'])
        end
    otherwise
        error(['Plotting side not specified correctly. ',...
            'Use sideToPlot = [''L'',''R'' or ''Both'']'])
end

dtptemp1 = dataToPlot1st;
dtptemp2 = dataToPlot2nd;
for j = 1:length(Stimuli) % for each OF stimulus
    for k = 1:length(neuronsToPlot)
        ni = neuronsToPlot(k);
        if ni ~= iH2
            switch EpochsFlipped(loc(j))
                case 1
                    dataToPlot1st{ni,j} = dtptemp1{ni,loc(j)};
                    dataToPlot2nd{ni,j} = dtptemp2{ni,loc(j)};
                case 2
                    dataToPlot2nd{ni,j} = dtptemp1{ni,loc(j)};
                    dataToPlot1st{ni,j} = dtptemp2{ni,loc(j)};
            end
        end
    end
end

frameNum = size(dataToPlot1st{1,1},1);
TimeAxisPD = linspace(0,frameNum/(minFrameRate/1000),frameNum);
iPeak = find(TimeAxisPD>4-peakbin/1000,1,'first'):...
    find(TimeAxisPD<4,1,'last');

% Calculate peak response, max normalize
Peaks1st = cell(size(dataToPlot1st));
Peaks2nd = cell(size(dataToPlot2nd));
Peaks1stNormM = cell(size(dataToPlot1st));
Peaks2ndNormM = cell(size(dataToPlot2nd));
PeaksDS = cell(size(dataToPlot1st));
PeaksDSNormM = cell(size(dataToPlot1st));
Maxes = nan(size(dataToPlot1st,1),1);
MaxesDS = cell(size(dataToPlot1st,1),1);
for n = 1:length(regions)
    my = zeros(1,length(Stimuli));
    myP = nan(1,size(dataToPlot1st{n,1},2));
    for ss = 1:length(stimToPlot)
        s = stimToPlot(ss);
        Peaks1st{n,s} = mean(dataToPlot1st{n,s}(iPeak,:));
        Peaks2nd{n,s} = mean(dataToPlot2nd{n,s}(iPeak,:));
        maxbl = max([nanmax(Peaks1st{n,s}),nanmax(Peaks2nd{n,s})]);
        if isempty(maxbl)
            maxbl = 0;
        end
        my(s) = max(my(s),maxbl);
        if contains(Stimuli{s},{'YawRightFlipped',...
                'Sideslip','Roll','Lift'})
            PeaksDS{n,s} = Peaks2nd{n,s}-Peaks1st{n,s};
        else
            if contains(Stimuli{s},{'Progressive',...
                    'ProgSide315','Yaw_','ProgSide45','Pitch'})
                PeaksDS{n,s} = Peaks1st{n,s}-Peaks2nd{n,s};
            end
        end
        maxbl = PeaksDS{n,s};
        myP = max(myP,maxbl);
    end
    Maxes(n) = max(my);
    MaxesDS{n} = myP;
    for ss = 1:length(stimToPlot)
        s = stimToPlot(ss);
        Peaks1stNormM{n,s} = Peaks1st{n,s}/Maxes(n);
        Peaks2ndNormM{n,s} = Peaks2nd{n,s}/Maxes(n);
        PeaksDSNormM{n,s} = PeaksDS{n,s}./MaxesDS{n};
    end
end

% plot the data
sp=0.5;
fs = 14; % font size for plot labels
alpha = 0.5; % transparency of the dots
ms = 32; % dot size of the mean data point

limsY = [-1 1.2];
f=figure('WindowState','maximized');

x = zeros(size(neuronsToPlot,2),size(stimToPlot,2));

my = zeros(1,size(neuronsToPlot,2));

for n = 1:length(neuronsToPlot)
    s = neuronsToPlot(n);
    iPlot = prepIndices{s} & gcIndices{s} & sideIndices{s};
    for k = 1:length(stimToPlot)
        m = stimToPlot(k);
        x(n,k) = (n-1)*length(stimToPlot) + 1 + (k-1)*sp;

        dataToPlot = PeaksDSNormM{s,m}(:,iPlot);
        scatter(ones(1,length(dataToPlot))*...
            x(n,k),dataToPlot,ms,...
            'MarkerFaceColor',cReg(s,:),'MarkerEdgeColor','none',...
            'MarkerFaceAlpha',alpha); hold on;
        plot(x(n,k),nanmean(dataToPlot),...
            'k.','MarkerSize',ms)
        eh = errorbar(x(n,k),nanmean(dataToPlot),...
            (nanstd(dataToPlot,0,2)/sqrt(sum(iPlot))));
        set(eh,'LineWidth',2,'Color','k')
        maxbl = nanmax(dataToPlot);
        if isempty(maxbl)
            maxbl = 0;
        end
        my(n) = max(my(n),maxbl);
        if k>1 && pairedLines
            for kk = 1:length(dataPrev)
                l = line([x(n,k-1),x(n,k)],[dataPrev(kk),...
                    dataToPlot(kk)]);
                set(l,'LineWidth',1,'Color','k','LineStyle','-')
            end
        end
        dataPrev = dataToPlot;
        dataToPlot(isnan(dataToPlot)) = [];
    end
end
limsX = xlim;
line(limsX,[0,0],'LineStyle','--','Color','k','LineWidth',1)
set(gca,'XTick',x(:,k)','XTickLabel',{regions{neuronsToPlot}},...
    'FontSize',14,'TickLabelInterpreter','none')
set(gca,'YLim',limsY,'YTick',linspace(-1,1,9),'FontSize',fs)
ylabel('DS Response (a.u.)')
set(gca,'XLim',limsX)

%% Figure 3e

% --------------------------- User Input ----------------------------------
% which neurons to plot. Defined by their index in the regions cell array.
% Plotted in the same order. Leave blank to plot all neurons.

neuronsToPlot = [8,4,11,5,12,3];

% Which calcium indicator(s) to plot. Defined by the genotype field in OF
% struct array. Leave blank to plot everything. Note that it is possible to
% pool data from multiple indicators (e.g. gcampsToPlot = {GC6f','GC7f'})
% Accepted arguments: {'GC6f','GC7f','sytGC7f'} Default: gcampsToPlot = {};
gcampsToPlot = {'sytGC7f','GC7f'};

% Which prep types to plot. Defined by the prepType field in OF struct
% array. Leave blank to plot all experiments.
% Accepted arguments: ['walking','immobilized'] Default: prepToPlot = '';
prepToPlot = 'immobilized';

% Neurons from which hemisphere to plot. Defined by the somaSide field in
% OF struct array. Type 'Both' to plot all neurons. In that case, responses
% from Right hemisphere will be flipped to match the Left hemisphere neural
% responses and plotted together.
% Accepted arguments: ['L','R','Both'] Default: sideToPlot = 'Both';
sideToPlot = 'Both';

% specify normalization type. Type 'Raw' for no normalization.
% 'Max' will normalize the data per ROI to the maximum response recorded in
% that ROI (see variable maxes).
% 'Zscore' will zscore the data per ROI based on the mean and standard
% deviation of the entire recording session done on that ROI (see variables
% mu and sigma). 'Znorm' normalizes Z-score data by the max Z-score
% observed per fly per neuron
% Accepted arguments: ['Zscore','Raw'] Default: 'Raw';
normType = 'Zscore';

% Specify plot type. Type 'Box' for plotting each fly with a different
% colored dot together with the population mean and sd as boxes.
% Currently only 'Box' is supported.
% Accepted arguments: ['Box']
plotType = 'Box';

% Specify grouping. Type 'neuron' for plotting all stimuli for each neuron
% grouped together. 'stimulus' will group all neurons for each stimulus
% Accepted arguments: ['neuron','stimulus']
groupType = 'neuron';

%----------------------------end of user input-----------------------------

peakbin = 500; % time bin for the averaged peak response (ms)
H2Flip = 0;
iH2 = find(contains(regions,'H2axon')); % index of H2
% index for recordings that will come from the right hemisphere
[~,loc] = ismember(StimuliFlipped,Stimuli);

clearvars l lh h stats
stm = {'Yaw_','Progressive','YawRightFlipped','RollInv','Pitch'};
auci = find(contains(Stimuli,stm) & BinocularIndices);
yi = auci(contains(Stimuli(auci),'Yaw_'));
yfi = auci(contains(Stimuli(auci),'YawRightFlipped'));
pri = auci(contains(Stimuli(auci),'Progressive'));
ssi = find(contains(Stimuli,'Sideslip'));
ps45i = find(contains(Stimuli,'ProgSide45'));

% Define the subset of data that will be plotted. This will be used in the
% plotSelectAUC function to plot only the specified subset of data as well
% as to calculate and display number of flies & ROIs being plotted

if isempty(neuronsToPlot)
    neuronsToPlot = 1:length(regions);
end

prepIndices = cell(length(regions),1);
gcIndices = cell(length(regions),1);
sideIndices = cell(length(regions),1);
prepIndicesFly = cell(length(regions),1);
gcIndicesFly = cell(length(regions),1);
sideIndicesFly = cell(length(regions),1);

if ~isempty(prepToPlot)
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        if isempty(prepIndices{r})
            prepIndices{r} = ...
                contains(PrepTypePerROI{r},prepToPlot);
            prepIndicesFly{r} = ...
                contains(PrepType{r},prepToPlot);
        else
            prepIndices{r} = (prepIndices{r} | ...
                contains(PrepTypePerROI{r},prepToPlot));
            prepIndicesFly{r} = (prepIndicesFly{r} | ...
                contains(PrepType{r},prepToPlot));
        end
    end
else
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        prepIndices{r} = true(size(PrepTypePerROI{r}));
        prepIndicesFly{r} = true(size(PrepType{r}));
    end
end
if ~isempty(gcampsToPlot)
    for p = 1:length(gcampsToPlot) %for each selected GCaMP type
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            if isempty(gcIndices{r})
                gcIndices{r} = ...
                    ismember(GCaMPTypePerROI{r},gcampsToPlot{p});
                gcIndicesFly{r} = ...
                    ismember(GCaMPType{r},gcampsToPlot{p});
            else
                gcIndices{r} = (gcIndices{r} | ...
                    ismember(GCaMPTypePerROI{r},gcampsToPlot{p}));
                gcIndicesFly{r} = (gcIndicesFly{r} | ...
                    ismember(GCaMPType{r},gcampsToPlot{p}));
            end
        end
    end
else
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        gcIndices{r} = true(size(GCaMPTypePerROI{r}));
        gcIndicesFly{r} = true(size(GCaMPType{r}));
    end
end
switch sideToPlot
    case 'Both'
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            sideIndices{r} = true(size(SomaSidePerROI{r}));
            sideIndicesFly{r} = true(size(SomaSide{r}));
        end
        switch normType
            case 'Zscore'
                dataToPlot1st = OFRes_1st_Z;
                dataToPlot2nd = OFRes_2nd_Z;
            case 'Max'
                dataToPlot1st = Auc_1st_normM;
                dataToPlot2nd = Auc_2nd_normM;
            case 'Raw'
                dataToPlot1st = Auc_1st_bs;
                dataToPlot2nd = Auc_2nd_bs;
            otherwise
                error(['Normalization type not specified correctly. ',...
                    'Use normType=[''Max'',''Zscore'', ',...
                    'or ''Raw'']'])
        end
    case {'L','R'}
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            if isempty(sideIndices{r})
                sideIndices{r} = ...
                    contains(SomaSidePerROI{r},sideToPlot);
                sideIndicesFly{r} = ...
                    contains(SomaSide{r},sideToPlot);
            else
                sideIndices{r} = (sideIndices{r} | ...
                    contains(SomaSidePerROI{r},sideToPlot));
                sideIndicesFly{r} = (sideIndicesFly{r} | ...
                    contains(SomaSide{r},sideToPlot));
            end
        end
        switch normType
            case 'Zscore'
                dataToPlot1st = Auc_1st_bs_Z;
                dataToPlot2nd = Auc_2nd_bs_Z;
            case 'Max'
                dataToPlot1st = Auc_1st_normM;
                dataToPlot2nd = Auc_2nd_normM;
            case 'Raw'
                dataToPlot1st = Auc_1st_bs;
                dataToPlot2nd = Auc_2nd_bs;
            otherwise
                error(['Normalization type not specified correctly. ',...
                    'Use normType=[''Max'',''Zscore'', ',...
                    'or ''Raw'']'])
        end
    otherwise
        error(['Plotting side not specified correctly. ',...
            'Use sideToPlot = [''L'',''R'' or ''Both'']'])
end

frameNum = size(dataToPlot1st{1,1},1);
TimeAxisPD = linspace(0,frameNum/(minFrameRate/1000),frameNum);
iPeak = find(TimeAxisPD>4-peakbin/1000,1,'first'):...
    find(TimeAxisPD<4,1,'last');

compNames = {'YawFlipped-Yaw','Prog-Yaw'};
ncomps = length(compNames);
sp=0.5;
x = zeros(ncomps,size(neuronsToPlot,2));
trIndex = cell(ncomps,size(neuronsToPlot,2));
f=figure('WindowState','maximized');
limsY = [-1.2 1.7];
my = zeros(1,ncomps);

for n = 1:ncomps
    for k = neuronsToPlot
        m = find(neuronsToPlot==k);
        x(n,m) = (n-1)*length(neuronsToPlot) + 1 + (m-1)*sp;
        iPlot = prepIndices{k} & gcIndices{k} & sideIndices{k};
        switch compNames{n}
            case 'Prog-Yaw'
                switch regions{k}
                    case {'H2axon'}
                        PKprog = mean(dataToPlot2nd{k,pri}(iPeak,iPlot))-mean(dataToPlot1st{k,pri}(iPeak,iPlot));
                        tPKprog = abs(mean(dataToPlot2nd{k,pri}(iPeak,iPlot))-mean(dataToPlot1st{k,pri}(iPeak,iPlot)));
                        PKyaw = mean(dataToPlot1st{k,yi}(iPeak,iPlot))-mean(dataToPlot2nd{k,yi}(iPeak,iPlot));
                        tPKyaw = abs(mean(dataToPlot1st{k,yi}(iPeak,iPlot))-mean(dataToPlot2nd{k,yi}(iPeak,iPlot)));
                        trIndex{n,m}=(PKprog-PKyaw)./(tPKprog+tPKyaw);
                    otherwise
                        PKprog = mean(dataToPlot1st{k,pri}(iPeak,iPlot))-mean(dataToPlot2nd{k,pri}(iPeak,iPlot));
                        tPKprog = abs(mean(dataToPlot1st{k,pri}(iPeak,iPlot))-mean(dataToPlot2nd{k,pri}(iPeak,iPlot)));
                        PKyaw = mean(dataToPlot2nd{k,yi}(iPeak,iPlot))-mean(dataToPlot1st{k,yi}(iPeak,iPlot));
                        tPKyaw = abs(mean(dataToPlot2nd{k,yi}(iPeak,iPlot))-mean(dataToPlot1st{k,yi}(iPeak,iPlot)));
                        trIndex{n,m}=(PKprog-PKyaw)./(tPKprog+tPKyaw);
                end
            case 'YawFlipped-Yaw'
                switch regions{k}
                    case {'H2axon'}
                        PKyf = mean(dataToPlot1st{k,yfi}(iPeak,iPlot))-mean(dataToPlot2nd{k,yfi}(iPeak,iPlot));
                        tPKyf = abs(mean(dataToPlot2nd{k,yfi}(iPeak,iPlot))-mean(dataToPlot1st{k,yfi}(iPeak,iPlot)));
                        PKyaw = mean(dataToPlot1st{k,yi}(iPeak,iPlot))-mean(dataToPlot2nd{k,yi}(iPeak,iPlot));
                        tPKyaw = abs(mean(dataToPlot1st{k,yi}(iPeak,iPlot))-mean(dataToPlot2nd{k,yi}(iPeak,iPlot)));
                        trIndex{n,m}=(PKyf-PKyaw)./(tPKyf+tPKyaw);
                    case {'H2rnaxon'}
                        PKyf = mean(dataToPlot1st{k,yfi}(iPeak,iPlot))-mean(dataToPlot2nd{k,yfi}(iPeak,iPlot));
                        tPKyf = abs(mean(dataToPlot2nd{k,yfi}(iPeak,iPlot))-mean(dataToPlot1st{k,yfi}(iPeak,iPlot)));
                        PKyaw = mean(dataToPlot2nd{k,yi}(iPeak,iPlot))-mean(dataToPlot1st{k,yi}(iPeak,iPlot));
                        tPKyaw = abs(mean(dataToPlot1st{k,yi}(iPeak,iPlot))-mean(dataToPlot2nd{k,yi}(iPeak,iPlot)));
                        trIndex{n,m}=(PKyf-PKyaw)./(tPKyf+tPKyaw);
                    otherwise
                        PKyf = mean(dataToPlot2nd{k,yfi}(iPeak,iPlot))-mean(dataToPlot1st{k,yfi}(iPeak,iPlot));
                        tPKyf = abs(mean(dataToPlot1st{k,yfi}(iPeak,iPlot))-mean(dataToPlot2nd{k,yfi}(iPeak,iPlot)));
                        PKyaw = mean(dataToPlot2nd{k,yi}(iPeak,iPlot))-mean(dataToPlot1st{k,yi}(iPeak,iPlot));
                        tPKyaw = abs(mean(dataToPlot2nd{k,yi}(iPeak,iPlot))-mean(dataToPlot1st{k,yi}(iPeak,iPlot)));
                        trIndex{n,m}=(PKyf-PKyaw)./(tPKyf+tPKyaw);
                end

            case 'Prog-Sideslip'
                switch regions{k}
                    case {'H2axon'}
                        PKss = mean(dataToPlot2nd{k,ssi}(iPeak,iPlot))-mean(dataToPlot1st{k,ssi}(iPeak,iPlot));
                        tPKss = abs(mean(dataToPlot2nd{k,ssi}(iPeak,iPlot))-mean(dataToPlot1st{k,ssi}(iPeak,iPlot)));
                        PKprog = mean(dataToPlot2nd{k,pri}(iPeak,iPlot))-mean(dataToPlot1st{k,pri}(iPeak,iPlot));
                        tPKprog = abs(mean(dataToPlot2nd{k,pri}(iPeak,iPlot))-mean(dataToPlot1st{k,pri}(iPeak,iPlot)));
                        trIndex{n,m}=(PKprog-PKss)./(tPKprog+tPKss);
                    otherwise
                        PKss = mean(dataToPlot1st{k,ssi}(iPeak,iPlot))-mean(dataToPlot2nd{k,ssi}(iPeak,iPlot));
                        tPKss = abs(mean(dataToPlot1st{k,ssi}(iPeak,iPlot))-mean(dataToPlot2nd{k,ssi}(iPeak,iPlot)));
                        PKprog = mean(dataToPlot1st{k,pri}(iPeak,iPlot))-mean(dataToPlot2nd{k,pri}(iPeak,iPlot));
                        tPKprog = abs(mean(dataToPlot1st{k,pri}(iPeak,iPlot))-mean(dataToPlot2nd{k,pri}(iPeak,iPlot)));
                        trIndex{n,m}=(PKprog-PKss)./(tPKprog+tPKss);
                end

            case 'ProgSide45-Sideslip'
                switch regions{k}
                    case {'H2axon'}
                        PKprog45 = mean(dataToPlot2nd{k,ps45i}(iPeak,iPlot))-mean(dataToPlot1st{k,ps45i}(iPeak,iPlot));
                        tPKprog45 = abs(mean(dataToPlot2nd{k,ps45i}(iPeak,iPlot))-mean(dataToPlot1st{k,ps45i}(iPeak,iPlot)));
                        PKss = mean(dataToPlot2nd{k,ssi}(iPeak,iPlot))-mean(dataToPlot1st{k,ssi}(iPeak,iPlot));
                        tPKss = abs(mean(dataToPlot2nd{k,ssi}(iPeak,iPlot))-mean(dataToPlot1st{k,ssi}(iPeak,iPlot)));
                        trIndex{n,m}=(PKprog45-PKss)./(tPKprog45+tPKss);
                    case {'H2rnaxon'}
                        PKprog45 = mean(dataToPlot1st{k,ps45i}(iPeak,iPlot))-mean(dataToPlot2nd{k,ps45i}(iPeak,iPlot));
                        tPKprog45 = abs(mean(dataToPlot1st{k,ps45i}(iPeak,iPlot))-mean(dataToPlot2nd{k,ps45i}(iPeak,iPlot)));
                        PKss = mean(dataToPlot1st{k,ssi}(iPeak,iPlot))-mean(dataToPlot2nd{k,ssi}(iPeak,iPlot));
                        tPKss = abs(mean(dataToPlot1st{k,ssi}(iPeak,iPlot))-mean(dataToPlot2nd{k,ssi}(iPeak,iPlot)));
                        trIndex{n,m}=(PKprog45-PKss)./(tPKprog45+tPKss);
                    otherwise
                        PKprog45 = mean(dataToPlot1st{k,ps45i}(iPeak,iPlot))-mean(dataToPlot2nd{k,ps45i}(iPeak,iPlot));
                        tPKprog45 = abs(mean(dataToPlot1st{k,ps45i}(iPeak,iPlot))-mean(dataToPlot2nd{k,ps45i}(iPeak,iPlot)));
                        PKss = mean(dataToPlot1st{k,ssi}(iPeak,iPlot))-mean(dataToPlot2nd{k,ssi}(iPeak,iPlot));
                        tPKss = abs(mean(dataToPlot1st{k,ssi}(iPeak,iPlot))-mean(dataToPlot2nd{k,ssi}(iPeak,iPlot)));
                        trIndex{n,m}=(PKprog45-PKss)./(tPKprog45+tPKss);
                end
        end
        h(n,m) = notBoxPlot(trIndex{n,m},x(n,m));hold on
        set(h(n,m).data,'MarkerFaceColor',cReg(k,:))
        % save the maximum y values for significance plot alignment
        maxbl = nanmax(trIndex{n,m});
        if isempty(maxbl)
            maxbl = 0;
        end
        my(n) = max(my(n),maxbl);
        trIndex{n,m}(isnan(trIndex{n,m})) = [];
    end
    % Perform Mann-Whitney-Wilcoxon non parametric test and plot
    % significance indicators if p<0.05
    nck = nchoosek(1:m,2);
    %This is all done just to help with the plots
    nck = [nck; nck(end,:)];
    nph = size(nck,1) - 1; % number of comparisons made (used for adjusting the
    % minimum P value required for significance - a.k.a. Bonferroni correction)
    sameComparison = 1;
    for ii = 1 : size(nck,1) - 1
        if nck(ii+1) == nck(ii)
            sameComparison = sameComparison + 1;
        else
            sameComparison = 1;
        end
    end
    for ii = 1 : size(nck,1) - 1
        stats{ii} = mwwtest(trIndex{n,nck(ii,1)},trIndex{n,nck(ii,2)},0);
        % disp([compNames{n},'-',regions{neuronsToPlot(nck(ii,1))},'-',...
        %     regions{neuronsToPlot(nck(ii,2))},'-',num2str(stats{ii}.p(2))])
        if stats{ii}.p < 0.05/nph
            % uncomment below if you want to display the p-values per pair
            disp([compNames{n},'-',regions{neuronsToPlot(nck(ii,1))},'-',...
                regions{neuronsToPlot(nck(ii,2))},'-',num2str(stats{ii}.p(2))])
            sigline([x(n,nck(ii,1)),x(n,nck(ii,2))],[],my(n)+my(n)*0.05*ii);
        end
    end
end
limsX = xlim;
line(limsX,[0,0],'LineStyle','--','Color','k','LineWidth',1)
set(gca,'XTick',x(:,end)','XTickLabel',compNames,'FontSize',18,...
    'YLim',limsY)
set(gca,'FontSize',14,'YLim',limsY)
ylabel('Discrimination Index (peak)')

for m=1:length(neuronsToPlot)
    hold on
    l(m)=plot(nan,nan,'Color',cReg(neuronsToPlot(m),:),'LineWidth',4);
end
set(gca,'XLim',limsX)
lh=legend(l,regions(neuronsToPlot),'Box','off','FontSize',24,'Location',...
    'Best');
clearvars l lh h

%% Figure 4c

% --------------------------- User Input ----------------------------------
% which OF stimuli to plot. Defined by their index in the Stimuli cell
% array. Plotted in the same order. Leave blank to plot all stimuli.

stimToPlot = [15 14 17 16]; % Fig 4c patterns

% which neurons to plot. Defined by their index in the regions cell array.
% Plotted in the same order. Leave blank to plot all neurons.

neuronsToPlot = [9,10]; % Fig 4c neurons

% Which calcium indicator(s) to plot. Defined by the genotype field in OF
% struct array. Leave blank to plot everything. Note that it is possible to
% pool data from multiple indicators (e.g. gcampsToPlot = {GC6f','GC7f'})
% Accepted arguments: {'GC6f','GC7f','sytGC7f','GC8f'}
% Default: gcampsToPlot = {};
gcampsToPlot = {'GC6f'};

% Which prep types to plot. Defined by the prepType field in OF struct
% array. Leave blank to plot all experiments. Keep in mind that all of the
% OF sensitivity experiments in walking prep were done when the fly was
% stationary and not moving/grooming/struggling etc.
% Accepted arguments: ['walking','immobilized'] Default: prepToPlot = '';
prepToPlot = 'immobilized';

% Neurons from which hemisphere to plot. Defined by the somaSide field in
% OF struct array. Type 'Both' to plot all neurons. In that case, responses
% from Right hemisphere will be flipped to match the Left hemisphere neural
% responses and plotted together.
% Accepted arguments: ['L','R','Both'] Default: sideToPlot = 'Both';
sideToPlot = 'Both';

% specify normalization type. Type 'DFF' for (F-Fo)/Fo normalization.
% 'Max' will normalize the data per ROI to the maximum DFF response
% recorded in that ROI (see variable maxes).
% 'Zscore' will z-score the data per ROI based on the baseline mean and
% standard deviation for each trial done on that ROI.
% 'Znorm' normalizes Z-score data by the max Z-score recorded from each ROI
% Accepted arguments: ['DFF','Max','Zscore','Znorm']
normType = 'Zscore';

% Specify plot type. Type 'Line' for plotting each fly with a different
% colored line together with the mean trace. 'Shaded' will plot the std as
% a shaded area around the mean trace instead.
% Accepted arguments: ['Line','Shaded']
plotType = 'Shaded';

% Define Y-axis limits for each normalization type.
switch normType
    case 'Max'
        limsY = [-0.5 1.1]; % Y-axis limits
    case 'DFF'
        limsY = [-3 10]; % Y-axis limits
    case 'Zscore'
        switch plotType
            case 'Line'
                limsY = [-6 50]; % Y-axis limits
            case 'Shaded'
                limsY = [-2 38]; % Y-axis limits
        end
    case 'Znorm'
        limsY = [-0.4 1.05]; % Y-axis limits
end
%----------------------------end of user input-----------------------------

H2flip = 1;
H2i = find(contains(regions,'H2axon')); % index of H2

% Define the subset of data that will be plotted. This will be used in the
% plotSelectDFF function to plot only the specified subset of data as well
% as to calculate and display number of flies & ROIs being plotted
if isempty(stimToPlot)
    stimToPlot = 1:length(Stimuli);
end
if isempty(neuronsToPlot)
    neuronsToPlot = 1:length(regions);
end

prepIndices = cell(length(regions),1);
gcIndices = cell(length(regions),1);
sideIndices = cell(length(regions),1);
prepIndicesFly = cell(length(regions),1);
gcIndicesFly = cell(length(regions),1);
sideIndicesFly = cell(length(regions),1);

if ~isempty(prepToPlot)
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        if isempty(prepIndices{r})
            prepIndices{r} = ...
                contains(PrepTypePerROI{r},prepToPlot);
            prepIndicesFly{r} = ...
                contains(PrepType{r},prepToPlot);
        else
            prepIndices{r} = (prepIndices{r} | ...
                contains(PrepTypePerROI{r},prepToPlot));
            prepIndicesFly{r} = (prepIndicesFly{r} | ...
                contains(PrepType{r},prepToPlot));
        end
    end
else
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        prepIndices{r} = true(size(PrepTypePerROI{r}));
        prepIndicesFly{r} = true(size(PrepType{r}));
    end
end
if ~isempty(gcampsToPlot)
    for p = 1:length(gcampsToPlot) %for each selected GCaMP type
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            if isempty(gcIndices{r})
                gcIndices{r} = ...
                    ismember(GCaMPTypePerROI{r},gcampsToPlot{p});
                gcIndicesFly{r} = ...
                    ismember(GCaMPType{r},gcampsToPlot{p});
            else
                gcIndices{r} = (gcIndices{r} | ...
                    ismember(GCaMPTypePerROI{r},gcampsToPlot{p}));
                gcIndicesFly{r} = (gcIndicesFly{r} | ...
                    ismember(GCaMPType{r},gcampsToPlot{p}));
            end
        end
    end
else
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        gcIndices{r} = true(size(GCaMPTypePerROI{r}));
        gcIndicesFly{r} = true(size(GCaMPType{r}));
    end
end

switch sideToPlot
    case 'Both'
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            sideIndices{r} = true(size(SomaSidePerROI{r}));
            sideIndicesFly{r} = true(size(SomaSide{r}));
        end
        switch normType
            case 'Zscore'
                dataToPlot1st = OFRes_1st_Z_R;
                dataToPlot2nd = OFRes_2nd_Z_R;
                if H2flip
                    dataToPlot1st(H2i,:) = OFRes_1st_Z(H2i,:);
                    dataToPlot2nd(H2i,:) = OFRes_2nd_Z(H2i,:);
                end
            case 'Max'
                dataToPlot1st = OFRes_1st_normM;
                dataToPlot2nd = OFRes_2nd_normM;
            case 'DFF'
                dataToPlot1st = OFRes_1st;
                dataToPlot2nd = OFRes_2nd;
            case 'Znorm'
                dataToPlot1st = OFRes_1st_Znorm;
                dataToPlot2nd = OFRes_2nd_Znorm;
            otherwise
                error(['Normalization type not specified correctly. ',...
                    'Use normType=[''Max'',''Zscore'',''Znorm'' ',...
                    'or ''Raw'']'])
        end
    case {'L','R'}
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            if isempty(sideIndices{r})
                sideIndices{r} = ...
                    contains(SomaSidePerROI{r},sideToPlot);
                sideIndicesFly{r} = ...
                    contains(SomaSide{r},sideToPlot);
            else
                sideIndices{r} = (sideIndices{r} | ...
                    contains(SomaSidePerROI{r},sideToPlot));
                sideIndicesFly{r} = (sideIndicesFly{r} | ...
                    contains(SomaSide{r},sideToPlot));
            end
        end
        switch normType
            case 'Zscore'
                dataToPlot1st = OFRes_Raw_1st_Z;
                dataToPlot2nd = OFRes_Raw_2nd_Z;
            case 'Max'
                dataToPlot1st = OFRes_Raw_1st_normM;
                dataToPlot2nd = OFRes_Raw_2nd_normM;
            case 'Raw'
                dataToPlot1st = OFRes_Raw_1st;
                dataToPlot2nd = OFRes_Raw_2nd;
            case 'Znorm'
                dataToPlot1st = OFRes_1st_Znorm;
                dataToPlot2nd = OFRes_2nd_Znorm;
            otherwise
                error(['Normalization type not specified correctly. ',...
                    'Use normType=[''Max'',''Zscore'',''Znorm'' ',...
                    'or ''Raw'']'])
        end
    otherwise
        error(['Plotting side not specified correctly. ',...
            'Use sideToPlot = [''L'',''R'' or ''Both'']'])
end

plotSelectDFF(dataToPlot1st,dataToPlot2nd,neuronsToPlot,...
    stimToPlot,minFrameRate,cReg,prepIndices,prepIndicesFly,...
    gcIndices,gcIndicesFly,sideIndices,sideIndicesFly,Stimuli,...
    regions,limsY,plotType,normType);

%% Figure 4d
% --------------------------- User Input ----------------------------------
neuronsToPlot = [9,10];
% Specify plot type. Type 'Line' for plotting each fly with a different
% colored line together with the mean trace. 'Shaded' will plot the std as
% a shaded area around the mean trace instead.
% Accepted arguments: ['Line','Shaded']
plotType = 'Shaded';
%----------------------------end of user input-----------------------------

% Translational stimuli used (as of 210421)
StimTR = {'Sideslip','ProgSide315','Progressive','ProgSide45'};

% Find the corresponding indices for the stimuli defined in StimTR in the
% Stimuli cell array
iTR = nan(1,length(StimTR));
for i = 1:length(iTR)
    iTR(i)=find(contains(Stimuli,StimTR(i))&~contains(Stimuli,'Left') & ...
        ~contains(Stimuli,'Right'));
end

lmY = [-3.5 23]; % Y-axis limits for bIPS-RdlWT
cRect = [0 0.5 0 0.1]; % color of the rectangle depicting stimulus motion
f=figure('WindowState','maximized'); hold on
for i = 1:length(neuronsToPlot) % for each unique neuron-region combination
    n = neuronsToPlot(i);
    frameNum = size(OFRes_1st_Z_R{n,1},1);
    TimeAxisPD = linspace(0,frameNum/(minFrameRate/1000),frameNum);
    %     limsX = [min(TimeAxisPD) ceil(max(TimeAxisPD))];
    limsX = [min(TimeAxisPD)+1 ceil(max(TimeAxisPD))-0.5];

    for j = 1:length(StimTR) % for each OF stimulus
        if ~isempty(OFRes_1st_Z_R{n,iTR(j)}) % if there is imaging data
            CounterTRn = sum(~isnan(OFRes_1st_Z_R{n,iTR(j)}(1,:)));
            switch StimTR{j}
                case 'Sideslip'
                    subplot(3,3,6)
                case 'ProgSide315'
                    subplot(3,3,1)
                case 'Progressive'
                    subplot(3,3,2)
                case 'ProgSide45'
                    subplot(3,3,3)
            end
            % plot individual neurites, and mean across neurites, 1st epoch
            switch plotType
                case 'Line'
                    plot(TimeAxisPD,OFRes_1st_Z_R{n,iTR(j)}); hold on;
                    plot(TimeAxisPD,nanmean(OFRes_1st_Z_R{n,iTR(j)},2),...
                        'Color',cReg(n,:),'LineWidth',4);
                case 'Shaded'
                    h=shadedErrorBar(TimeAxisPD,...
                        nanmean(OFRes_1st_Z_R{n,iTR(j)},2),...
                        nanstd(OFRes_1st_Z_R{n,iTR(j)},0,2)/...
                        sqrt(CounterTRn));hold on
                    set(h.mainLine,'Color',cReg(n,:),...
                        'LineWidth',4)
                    set(h.patch,'FaceColor',cReg(n,:),...
                        'FaceAlpha',0.2)
                    set(h.edge(1),'Color',cReg(n,:))
                    set(h.edge(2),'Color',cReg(n,:))
                    clearvars h
            end

            rectangle('Position',[2 lmY(1) 2 lmY(2)-lmY(1)],...
                'FaceColor', cRect);
            set(gca,'YLim',lmY,'XLim',limsX)
            switch StimTR{j}
                case 'ProgSide315'
                    ylabel('Response (z-score)')
                    xlabel('Time (sec)')
                    set(gca,'FontSize',18)
                otherwise
                    axis off
            end
            line(limsX,[0 0],'Color','k','LineWidth',0.5,'LineStyle','--')
        end

        if  ~isempty(OFRes_2nd_Z_R{n,iTR(j)}) % if there is imaging data
            switch StimTR{j}
                case 'Sideslip'
                    subplot(3,3,4)
                case 'ProgSide315'
                    subplot(3,3,9)
                case 'Progressive'
                    subplot(3,3,8)
                case 'ProgSide45'
                    subplot(3,3,7)
            end
            % plot individual neurites, and mean across neurites, 2nd epoch
            switch plotType
                case 'Line'
                    plot(TimeAxisPD,OFRes_2nd_Z_R{n,iTR(j)}); hold on;
                    plot(TimeAxisPD,nanmean(OFRes_2nd_Z_R{n,iTR(j)},2),...
                        'color',cReg(n,:),'LineWidth',4);
                case 'Shaded'
                    h=shadedErrorBar(TimeAxisPD,...
                        nanmean(OFRes_2nd_Z_R{n,iTR(j)},2),...
                        nanstd(OFRes_2nd_Z_R{n,iTR(j)},0,2)/...
                        sqrt(CounterTRn));hold on
                    set(h.mainLine,'Color',cReg(n,:),...
                        'LineWidth',4)
                    set(h.patch,'FaceColor',cReg(n,:),...
                        'FaceAlpha',0.2)
                    set(h.edge(1),'Color',cReg(n,:))
                    set(h.edge(2),'Color',cReg(n,:))
                    clearvars h
            end
            rectangle('Position',[2 lmY(1) 2 lmY(2)-lmY(1)],...
                'FaceColor', cRect);
            set(gca,'YLim',lmY,'XLim',limsX)
            axis off
            line(limsX,[0 0],'Color','k','LineWidth',0.5,'LineStyle','--')
        end
    end
end

%% Figure 4e

% Here you can specify which calcium imaging responses (peak) you'd like to
% plot for visualization. You can select which regions, calcium indicators,
% soma side of the imaged neurons, prep type and OF stimuli (non ST).
% Finally, you can specify if and which normalization should be applied to
% the data

% --------------------------- User Input ----------------------------------
% which OF stimuli to plot. Defined by their index in the Stimuli cell
% array. Plotted in the same order. Leave blank to plot all stimuli.
stimToPlot = [14 15 17 16 13 7]; % Yaw YawFlipped Prog Sideslip

% which neurons to plot. Defined by their index in the regions cell array.
% Plotted in the same order. Leave blank to plot all neurons.
neuronsToPlot = [9,10];

% Which calcium indicator(s) to plot. Defined by the genotype field in OF
% struct array. Leave blank to plot everything. Note that it is possible to
% pool data from multiple indicators (e.g. gcampsToPlot = {GC6f','GC7f'})
% Accepted arguments: {'GC6f','GC7f','sytGC7f'} Default: gcampsToPlot = {};
gcampsToPlot = {'GC6f'};

% Which prep types to plot. Defined by the prepType field in OF struct
% array. Leave blank to plot all experiments.
% Accepted arguments: ['walking','immobilized'] Default: prepToPlot = '';
prepToPlot = 'immobilized';

% Neurons from which hemisphere to plot. Defined by the somaSide field in
% OF struct array. Type 'Both' to plot all neurons. In that case, responses
% from Right hemisphere will be flipped to match the Left hemisphere neural
% responses and plotted together.
% Accepted arguments: ['L','R','Both'] Default: sideToPlot = 'Both';
sideToPlot = 'Both';

% specify normalization type. Type 'Raw' for no normalization.
% 'Max' will normalize the data per ROI to the maximum response recorded in
% that ROI (see variable maxes).
% 'Zscore' will zscore the data per ROI based on the mean and standard
% deviation of the entire recording session done on that ROI (see variables
% mu and sigma). 'Znorm' normalizes Z-score data by the max Z-score
% observed per fly per neuron
% Accepted arguments: ['Zscore','Raw'] Default: 'Raw';
normType = 'Zscore';

% Specify grouping. Type 'neuron' for plotting all stimuli for each neuron
% grouped together. 'stimulus' will group all neurons for each stimulus
% Accepted arguments: ['neuron','stimulus']
groupType = 'stimulus';


%----------------------------end of user input-----------------------------

peakbin = 500; % time bin for the averaged peak response (ms)
iH2 = find(contains(regions,'H2axon')); % index of H2
% index for recordings that will come from the right hemisphere
[~,loc] = ismember(StimuliFlipped,Stimuli);

clearvars l lh h stats
stm = {'Yaw_','Progressive','YawRightFlipped','RollInv','Pitch'};
auci = find(contains(Stimuli,stm) & BinocularIndices);
yi = auci(contains(Stimuli(auci),'Yaw_')); % index of the stimulus that
% will be used for normalization (currently is binocular yaw)

% Define the subset of data that will be plotted. This will be used in the
% plotSelectAUC function to plot only the specified subset of data as well
% as to calculate and display number of flies & ROIs being plotted
if isempty(stimToPlot)
    stimToPlot = 1:length(Stimuli);
end
if isempty(neuronsToPlot)
    neuronsToPlot = 1:length(regions);
end

prepIndices = cell(length(regions),1);
gcIndices = cell(length(regions),1);
sideIndices = cell(length(regions),1);
prepIndicesFly = cell(length(regions),1);
gcIndicesFly = cell(length(regions),1);
sideIndicesFly = cell(length(regions),1);

if ~isempty(prepToPlot)
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        if isempty(prepIndices{r})
            prepIndices{r} = ...
                contains(PrepTypePerROI{r},prepToPlot);
            prepIndicesFly{r} = ...
                contains(PrepType{r},prepToPlot);
        else
            prepIndices{r} = (prepIndices{r} | ...
                contains(PrepTypePerROI{r},prepToPlot));
            prepIndicesFly{r} = (prepIndicesFly{r} | ...
                contains(PrepType{r},prepToPlot));
        end
    end
else
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        prepIndices{r} = true(size(PrepTypePerROI{r}));
        prepIndicesFly{r} = true(size(PrepType{r}));
    end
end
if ~isempty(gcampsToPlot)
    for p = 1:length(gcampsToPlot) %for each selected GCaMP type
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            if isempty(gcIndices{r})
                gcIndices{r} = ...
                    ismember(GCaMPTypePerROI{r},gcampsToPlot{p});
                gcIndicesFly{r} = ...
                    ismember(GCaMPType{r},gcampsToPlot{p});
            else
                gcIndices{r} = (gcIndices{r} | ...
                    ismember(GCaMPTypePerROI{r},gcampsToPlot{p}));
                gcIndicesFly{r} = (gcIndicesFly{r} | ...
                    ismember(GCaMPType{r},gcampsToPlot{p}));
            end
        end
    end
else
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        gcIndices{r} = true(size(GCaMPTypePerROI{r}));
        gcIndicesFly{r} = true(size(GCaMPType{r}));
    end
end
switch sideToPlot
    case 'Both'
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            sideIndices{r} = true(size(SomaSidePerROI{r}));
            sideIndicesFly{r} = true(size(SomaSide{r}));
        end
        switch normType
            case 'Zscore'
                dataToPlot1st = OFRes_1st_Z;
                dataToPlot2nd = OFRes_2nd_Z;
            case 'Max'
                dataToPlot1st = Auc_1st_normM;
                dataToPlot2nd = Auc_2nd_normM;
            case 'Raw'
                dataToPlot1st = Auc_1st_bs;
                dataToPlot2nd = Auc_2nd_bs;
            otherwise
                error(['Normalization type not specified correctly. ',...
                    'Use normType=[''Max'',''Zscore'', ',...
                    'or ''Raw'']'])
        end
    case {'L','R'}
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            if isempty(sideIndices{r})
                sideIndices{r} = ...
                    contains(SomaSidePerROI{r},sideToPlot);
                sideIndicesFly{r} = ...
                    contains(SomaSide{r},sideToPlot);
            else
                sideIndices{r} = (sideIndices{r} | ...
                    contains(SomaSidePerROI{r},sideToPlot));
                sideIndicesFly{r} = (sideIndicesFly{r} | ...
                    contains(SomaSide{r},sideToPlot));
            end
        end
        switch normType
            case 'Zscore'
                dataToPlot1st = Auc_1st_bs_Z;
                dataToPlot2nd = Auc_2nd_bs_Z;
            case 'Max'
                dataToPlot1st = Auc_1st_normM;
                dataToPlot2nd = Auc_2nd_normM;
            case 'Raw'
                dataToPlot1st = Auc_1st_bs;
                dataToPlot2nd = Auc_2nd_bs;
            otherwise
                error(['Normalization type not specified correctly. ',...
                    'Use normType=[''Max'',''Zscore'', ',...
                    'or ''Raw'']'])
        end
    otherwise
        error(['Plotting side not specified correctly. ',...
            'Use sideToPlot = [''L'',''R'' or ''Both'']'])
end


dtptemp1 = dataToPlot1st;
dtptemp2 = dataToPlot2nd;
for j = 1:length(Stimuli) % for each OF stimulus
    for k = 1:length(neuronsToPlot)
        ni = neuronsToPlot(k);
        if ni ~= iH2
            switch EpochsFlipped(loc(j))
                case 1
                    dataToPlot1st{ni,j} = dtptemp1{ni,loc(j)};
                    dataToPlot2nd{ni,j} = dtptemp2{ni,loc(j)};
                case 2
                    dataToPlot2nd{ni,j} = dtptemp1{ni,loc(j)};
                    dataToPlot1st{ni,j} = dtptemp2{ni,loc(j)};
            end
        end
    end
end

frameNum = size(dataToPlot1st{1,1},1);
TimeAxisPD = linspace(0,frameNum/(minFrameRate/1000),frameNum);
iPeak = find(TimeAxisPD>4-peakbin/1000,1,'first'):...
    find(TimeAxisPD<4,1,'last');

sp=0.5;
ms = 38; % dot size of the mean data point
f=figure('WindowState','maximized');
switch groupType
    case 'stimulus'
        x = zeros(size(stimToPlot,2),size(neuronsToPlot,2));
        DSI = cell(length(stimToPlot),size(neuronsToPlot,2));
        DSR = cell(length(stimToPlot),size(neuronsToPlot,2));
        limsY = [-20 170];
        my = zeros(1,length(stimToPlot));
        for n = 1:length(stimToPlot)
            s = stimToPlot(n);
            for k = neuronsToPlot
                m = find(neuronsToPlot==k);
                x(n,m) = (n-1)*length(neuronsToPlot) + 1 + (m-1)*sp;
                iPlot = prepIndices{k} & gcIndices{k} & sideIndices{k};
                if contains(Stimuli{s},{'Yaw_','ProgSide45','Roll'})
                    switch regions{k}
                        case {'H2axon'}
                            rNPD = mean(dataToPlot2nd{k,s}(iPeak,iPlot));
                            rPD = mean(dataToPlot1st{k,s}(iPeak,iPlot));
                        otherwise
                            rPD = mean(dataToPlot1st{k,s}(iPeak,iPlot));
                            rNPD = mean(dataToPlot2nd{k,s}(iPeak,iPlot));
                    end
                else
                    if contains(Stimuli{stimToPlot(n)},{'Progressive','ProgSide315'})
                        switch regions{k}
                            case {'H2rnaxon','H2axon'}
                                rPD = mean(dataToPlot2nd{k,s}(iPeak,iPlot));
                                rNPD = mean(dataToPlot1st{k,s}(iPeak,iPlot));
                            otherwise
                                rPD = mean(dataToPlot1st{k,s}(iPeak,iPlot));
                                rNPD = mean(dataToPlot2nd{k,s}(iPeak,iPlot));
                        end
                    else
                        if contains(Stimuli{stimToPlot(n)},{'Pitch'})
                            switch regions{k}
                                case {'DNp15dendrite','DNp15-RdlFlpStop-WTdendrite',...
                                        'DNp15-RdlFlpStopdendrite'}
                                    rPD = mean(dataToPlot2nd{k,s}(iPeak,iPlot));
                                    rNPD = mean(dataToPlot1st{k,s}(iPeak,iPlot));
                                otherwise
                                    rPD = mean(dataToPlot1st{k,s}(iPeak,iPlot));
                                    rNPD = mean(dataToPlot2nd{k,s}(iPeak,iPlot));
                            end
                        else
                            if contains(Stimuli{stimToPlot(n)},{'Lift','YawRightFlipped','Sideslip'})
                                switch regions{k}
                                    case {'H2axon'}
                                        rNPD = mean(dataToPlot2nd{k,s}(iPeak,iPlot));
                                        rPD = mean(dataToPlot1st{k,s}(iPeak,iPlot));
                                    otherwise
                                        rNPD = mean(dataToPlot1st{k,s}(iPeak,iPlot));
                                        rPD = mean(dataToPlot2nd{k,s}(iPeak,iPlot));
                                end
                            end
                        end
                    end
                end
                DSI{n,m}=(rPD-rNPD)./(abs(rPD)+abs(rNPD));
                DSR{n,m}=(rPD-rNPD);
                h(n,m) = notBoxPlot(DSR{n,m},x(n,m));hold on
                set(h(n,m).data,'MarkerFaceColor',cReg(k,:))
                % save the maximum y values for significance plot alignment
                maxbl = nanmax(DSR{n,m});
                if isempty(maxbl)
                    maxbl = 0;
                end
                my(n) = max(my(n),maxbl);
                DSI{n,m}(isnan(DSI{n,m})) = [];
                DSR{n,m}(isnan(DSR{n,m})) = [];
            end
            % Perform Mann-Whitney-Wilcoxon non parametric test and plot
            % significance indicators if p<0.05
            nck = nchoosek(1:m,2);
            %This is all done just to help with the plots
            nck = [nck; nck(end,:)];
            nph = size(nck,1) - 1; % number of comparisons made (used for adjusting the
            % minimum P value required for significance - a.k.a. Bonferroni correction)
            sameComparison = 1;
            for ii = 1 : size(nck,1) - 1
                if nck(ii+1) == nck(ii)
                    sameComparison = sameComparison + 1;
                else
                    sameComparison = 1;
                end
            end
            for ii = 1 : size(nck,1) - 1
                stats{ii} = mwwtest(DSR{n,nck(ii,1)},DSR{n,nck(ii,2)},0);

                if stats{ii}.p < 0.05/nph
                    % uncomment below if you want to display the p-values per pair
                    disp(['n=',num2str(n),'-',regions{neuronsToPlot(nck(ii,1))},'-',...
                        regions{neuronsToPlot(nck(ii,2))},'-',num2str(stats{ii}.p(2))])
                    sigline([x(n,nck(ii,1)),x(n,nck(ii,2))],[],my(n)+my(n)*0.03*ii);
                end
            end
        end
        limsX = xlim;
        line(limsX,[0,0],'LineStyle','--','Color','k','LineWidth',1)
        set(gca,'XTick',x(:,2)','XTickLabel',{Stimuli{stimToPlot}},...
            'FontSize',14,'TickLabelInterpreter','none')
        ylabel('OFR (z-score)')
        
        for m=1:length(neuronsToPlot)
            hold on
            l(m)=plot(nan,nan,'Color',cReg(neuronsToPlot(m),:),'LineWidth',4);
        end
        set(gca,'XLim',limsX,'YLim',limsY)
        lh=legend(l,regions(neuronsToPlot),'Box','off','FontSize',24,'Location','Best');

    case 'neuron'
        x = zeros(size(neuronsToPlot,2),size(stimToPlot,2));
        DSI = cell(size(neuronsToPlot,2),size(stimToPlot,2));
        DSR = cell(size(neuronsToPlot,2),size(stimToPlot,2));
        limsY = [-20 170];
        my = zeros(1,size(neuronsToPlot,2));

        for n = 1:length(neuronsToPlot)
            s = neuronsToPlot(n);
            for k = 1:length(stimToPlot)
                m = stimToPlot(k);
                x(n,k) = (n-1)*length(stimToPlot) + 1 + (k-1)*sp;
                iPlot = prepIndices{s} & gcIndices{s} & sideIndices{s};
                if contains(Stimuli{m},{'YawRightFlipped','Roll','Sideslip','Lift'})
                    rNPD = mean(dataToPlot1st{s,m}(iPeak,iPlot));
                    rPD = mean(dataToPlot2nd{s,m}(iPeak,iPlot));
                else
                    if contains(Stimuli{m},{'Progressive','ProgSide45','ProgSide315','Yaw_','Pitch'})
                        rPD = mean(dataToPlot1st{s,m}(iPeak,iPlot));
                        rNPD = mean(dataToPlot2nd{s,m}(iPeak,iPlot));
                    end
                end
                DSI{n,k}=(rPD-rNPD)./(abs(rPD)+abs(rNPD));
                DSR{n,k}=(rPD-rNPD);
                h(n,k) = notBoxPlot(DSR{n,k},x(n,k));hold on
                set(h(n,k).data,'MarkerFaceColor',cReg(s,:))

                maxbl = nanmax(DSR{n,k});
                if isempty(maxbl)
                    maxbl = 0;
                end
                my(n) = max(my(n),maxbl);
                DSI{n,k}(isnan(DSI{n,k})) = [];
                DSR{n,k}(isnan(DSR{n,k})) = [];

            end
            % Perform Mann-Whitney-Wilcoxon non parametric test and plot
            % significance indicators if p<0.05
            nck = nchoosek(1:k,2);
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
                stats{ii} = mwwtest(DSR{n,nck(ii,1)},DSR{n,nck(ii,2)},0);
                if stats{ii}.p < 0.05
                    sigline([x(n,nck(ii,1)),x(n,nck(ii,2))],[],my(n)+my(n)*0.03*ii);
                    % uncomment below if you want to display the p-values per pair
                    disp(['n=',num2str(n),'-',regions{neuronsToPlot(n)},'-',...
                        Stimuli{stimToPlot(nck(ii,2))},'-',num2str(stats{ii}.p(2))])
                end
            end
        end
        limsX = xlim;
        line(limsX,[0,0],'LineStyle','--','Color','k','LineWidth',1)
        set(gca,'XTick',x(:,k)','XTickLabel',{regions{neuronsToPlot}},...
            'FontSize',14,'TickLabelInterpreter','none')
        set(gca,'YLim',limsY)
        ylabel('OFR (z-score)')
        for m=1:length(neuronsToPlot)
            hold on
            l(m)=plot(nan,nan,'Color',cReg(neuronsToPlot(n),:),'LineWidth',4);
        end
        set(gca,'XLim',limsX,'YLim',limsY)
        lh=legend(l,Stimuli(stimToPlot),'Box','off','FontSize',20,...
            'Location','Best','Interpreter','none');
    otherwise
        error(['Grouping type not specified correctly. ',...
            'Use groupType = [''neuron'' or ''stimulus'']'])
end

%% Figure 4f

% --------------------------- User Input ----------------------------------
% which neurons to plot. Defined by their index in the regions cell array.
% Plotted in the same order. Leave blank to plot all neurons.

neuronsToPlot = [9,10];

% Which calcium indicator(s) to plot. Defined by the genotype field in OF
% struct array. Leave blank to plot everything. Note that it is possible to
% pool data from multiple indicators (e.g. gcampsToPlot = {GC6f','GC7f'})
% Accepted arguments: {'GC6f','GC7f','sytGC7f'} Default: gcampsToPlot = {};
gcampsToPlot = {'GC6f'};

% Which prep types to plot. Defined by the prepType field in OF struct
% array. Leave blank to plot all experiments.
% Accepted arguments: ['walking','immobilized'] Default: prepToPlot = '';
prepToPlot = 'immobilized';

% Neurons from which hemisphere to plot. Defined by the somaSide field in
% OF struct array. Type 'Both' to plot all neurons. In that case, responses
% from Right hemisphere will be flipped to match the Left hemisphere neural
% responses and plotted together.
% Accepted arguments: ['L','R','Both'] Default: sideToPlot = 'Both';
sideToPlot = 'Both';

% specify normalization type. Type 'Raw' for no normalization.
% 'Max' will normalize the data per ROI to the maximum response recorded in
% that ROI (see variable maxes).
% 'Zscore' will zscore the data per ROI based on the mean and standard
% deviation of the entire recording session done on that ROI (see variables
% mu and sigma). 'Znorm' normalizes Z-score data by the max Z-score
% observed per fly per neuron
% Accepted arguments: ['Zscore','Raw'] Default: 'Raw';
normType = 'Zscore';

% Specify plot type. Type 'Box' for plotting each fly with a different
% colored dot together with the population mean and sd as boxes.
% Currently only 'Box' is supported.
% Accepted arguments: ['Box']
plotType = 'Box';

% Specify grouping. Type 'neuron' for plotting all stimuli for each neuron
% grouped together. 'stimulus' will group all neurons for each stimulus
% Accepted arguments: ['neuron','stimulus']
groupType = 'neuron';

%----------------------------end of user input-----------------------------

peakbin = 500; % time bin for the averaged peak response (ms)
H2Flip = 0;
iH2 = find(contains(regions,'H2axon')); % index of H2
% index for recordings that will come from the right hemisphere
[~,loc] = ismember(StimuliFlipped,Stimuli);

clearvars l lh h stats
stm = {'Yaw_','Progressive','YawRightFlipped','RollInv','Pitch'};
auci = find(contains(Stimuli,stm) & BinocularIndices);
yi = auci(contains(Stimuli(auci),'Yaw_'));
yfi = auci(contains(Stimuli(auci),'YawRightFlipped'));
pri = auci(contains(Stimuli(auci),'Progressive'));
ssi = find(contains(Stimuli,'Sideslip'));
ps45i = find(contains(Stimuli,'ProgSide45'));

% Define the subset of data that will be plotted. This will be used in the
% plotSelectAUC function to plot only the specified subset of data as well
% as to calculate and display number of flies & ROIs being plotted

if isempty(neuronsToPlot)
    neuronsToPlot = 1:length(regions);
end

prepIndices = cell(length(regions),1);
gcIndices = cell(length(regions),1);
sideIndices = cell(length(regions),1);
prepIndicesFly = cell(length(regions),1);
gcIndicesFly = cell(length(regions),1);
sideIndicesFly = cell(length(regions),1);

if ~isempty(prepToPlot)
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        if isempty(prepIndices{r})
            prepIndices{r} = ...
                contains(PrepTypePerROI{r},prepToPlot);
            prepIndicesFly{r} = ...
                contains(PrepType{r},prepToPlot);
        else
            prepIndices{r} = (prepIndices{r} | ...
                contains(PrepTypePerROI{r},prepToPlot));
            prepIndicesFly{r} = (prepIndicesFly{r} | ...
                contains(PrepType{r},prepToPlot));
        end
    end
else
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        prepIndices{r} = true(size(PrepTypePerROI{r}));
        prepIndicesFly{r} = true(size(PrepType{r}));
    end
end
if ~isempty(gcampsToPlot)
    for p = 1:length(gcampsToPlot) %for each selected GCaMP type
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            if isempty(gcIndices{r})
                gcIndices{r} = ...
                    ismember(GCaMPTypePerROI{r},gcampsToPlot{p});
                gcIndicesFly{r} = ...
                    ismember(GCaMPType{r},gcampsToPlot{p});
            else
                gcIndices{r} = (gcIndices{r} | ...
                    ismember(GCaMPTypePerROI{r},gcampsToPlot{p}));
                gcIndicesFly{r} = (gcIndicesFly{r} | ...
                    ismember(GCaMPType{r},gcampsToPlot{p}));
            end
        end
    end
else
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        gcIndices{r} = true(size(GCaMPTypePerROI{r}));
        gcIndicesFly{r} = true(size(GCaMPType{r}));
    end
end
switch sideToPlot
    case 'Both'
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            sideIndices{r} = true(size(SomaSidePerROI{r}));
            sideIndicesFly{r} = true(size(SomaSide{r}));
        end
        switch normType
            case 'Zscore'
                dataToPlot1st = OFRes_1st_Z;
                dataToPlot2nd = OFRes_2nd_Z;
            case 'Max'
                dataToPlot1st = Auc_1st_normM;
                dataToPlot2nd = Auc_2nd_normM;
            case 'Raw'
                dataToPlot1st = Auc_1st_bs;
                dataToPlot2nd = Auc_2nd_bs;
            otherwise
                error(['Normalization type not specified correctly. ',...
                    'Use normType=[''Max'',''Zscore'', ',...
                    'or ''Raw'']'])
        end
    case {'L','R'}
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            if isempty(sideIndices{r})
                sideIndices{r} = ...
                    contains(SomaSidePerROI{r},sideToPlot);
                sideIndicesFly{r} = ...
                    contains(SomaSide{r},sideToPlot);
            else
                sideIndices{r} = (sideIndices{r} | ...
                    contains(SomaSidePerROI{r},sideToPlot));
                sideIndicesFly{r} = (sideIndicesFly{r} | ...
                    contains(SomaSide{r},sideToPlot));
            end
        end
        switch normType
            case 'Zscore'
                dataToPlot1st = Auc_1st_bs_Z;
                dataToPlot2nd = Auc_2nd_bs_Z;
            case 'Max'
                dataToPlot1st = Auc_1st_normM;
                dataToPlot2nd = Auc_2nd_normM;
            case 'Raw'
                dataToPlot1st = Auc_1st_bs;
                dataToPlot2nd = Auc_2nd_bs;
            otherwise
                error(['Normalization type not specified correctly. ',...
                    'Use normType=[''Max'',''Zscore'', ',...
                    'or ''Raw'']'])
        end
    otherwise
        error(['Plotting side not specified correctly. ',...
            'Use sideToPlot = [''L'',''R'' or ''Both'']'])
end

frameNum = size(dataToPlot1st{1,1},1);
TimeAxisPD = linspace(0,frameNum/(minFrameRate/1000),frameNum);
iPeak = find(TimeAxisPD>4-peakbin/1000,1,'first'):...
    find(TimeAxisPD<4,1,'last');

compNames = {'YawFlipped-Yaw','Prog-Yaw'};
ncomps = length(compNames);
sp=0.5;
x = zeros(ncomps,size(neuronsToPlot,2));
trIndex = cell(ncomps,size(neuronsToPlot,2));
f=figure('WindowState','maximized');
limsY = [-1.2 1.7];
my = zeros(1,ncomps);

for n = 1:ncomps
    for k = neuronsToPlot
        m = find(neuronsToPlot==k);
        x(n,m) = (n-1)*length(neuronsToPlot) + 1 + (m-1)*sp;
        iPlot = prepIndices{k} & gcIndices{k} & sideIndices{k};
        switch compNames{n}
            case 'Prog-Yaw'
                switch regions{k}
                    case {'H2axon'}
                        PKprog = mean(dataToPlot2nd{k,pri}(iPeak,iPlot))-mean(dataToPlot1st{k,pri}(iPeak,iPlot));
                        tPKprog = abs(mean(dataToPlot2nd{k,pri}(iPeak,iPlot))-mean(dataToPlot1st{k,pri}(iPeak,iPlot)));
                        PKyaw = mean(dataToPlot1st{k,yi}(iPeak,iPlot))-mean(dataToPlot2nd{k,yi}(iPeak,iPlot));
                        tPKyaw = abs(mean(dataToPlot1st{k,yi}(iPeak,iPlot))-mean(dataToPlot2nd{k,yi}(iPeak,iPlot)));
                        trIndex{n,m}=(PKprog-PKyaw)./(tPKprog+tPKyaw);
                    otherwise
                        PKprog = mean(dataToPlot1st{k,pri}(iPeak,iPlot))-mean(dataToPlot2nd{k,pri}(iPeak,iPlot));
                        tPKprog = abs(mean(dataToPlot1st{k,pri}(iPeak,iPlot))-mean(dataToPlot2nd{k,pri}(iPeak,iPlot)));
                        PKyaw = mean(dataToPlot2nd{k,yi}(iPeak,iPlot))-mean(dataToPlot1st{k,yi}(iPeak,iPlot));
                        tPKyaw = abs(mean(dataToPlot2nd{k,yi}(iPeak,iPlot))-mean(dataToPlot1st{k,yi}(iPeak,iPlot)));
                        trIndex{n,m}=(PKprog-PKyaw)./(tPKprog+tPKyaw);
                end
            case 'YawFlipped-Yaw'
                switch regions{k}
                    case {'H2axon'}
                        PKyf = mean(dataToPlot1st{k,yfi}(iPeak,iPlot))-mean(dataToPlot2nd{k,yfi}(iPeak,iPlot));
                        tPKyf = abs(mean(dataToPlot2nd{k,yfi}(iPeak,iPlot))-mean(dataToPlot1st{k,yfi}(iPeak,iPlot)));
                        PKyaw = mean(dataToPlot1st{k,yi}(iPeak,iPlot))-mean(dataToPlot2nd{k,yi}(iPeak,iPlot));
                        tPKyaw = abs(mean(dataToPlot1st{k,yi}(iPeak,iPlot))-mean(dataToPlot2nd{k,yi}(iPeak,iPlot)));
                        trIndex{n,m}=(PKyf-PKyaw)./(tPKyf+tPKyaw);
                    case {'H2rnaxon'}
                        PKyf = mean(dataToPlot1st{k,yfi}(iPeak,iPlot))-mean(dataToPlot2nd{k,yfi}(iPeak,iPlot));
                        tPKyf = abs(mean(dataToPlot2nd{k,yfi}(iPeak,iPlot))-mean(dataToPlot1st{k,yfi}(iPeak,iPlot)));
                        PKyaw = mean(dataToPlot2nd{k,yi}(iPeak,iPlot))-mean(dataToPlot1st{k,yi}(iPeak,iPlot));
                        tPKyaw = abs(mean(dataToPlot1st{k,yi}(iPeak,iPlot))-mean(dataToPlot2nd{k,yi}(iPeak,iPlot)));
                        trIndex{n,m}=(PKyf-PKyaw)./(tPKyf+tPKyaw);
                    otherwise
                        PKyf = mean(dataToPlot2nd{k,yfi}(iPeak,iPlot))-mean(dataToPlot1st{k,yfi}(iPeak,iPlot));
                        tPKyf = abs(mean(dataToPlot1st{k,yfi}(iPeak,iPlot))-mean(dataToPlot2nd{k,yfi}(iPeak,iPlot)));
                        PKyaw = mean(dataToPlot2nd{k,yi}(iPeak,iPlot))-mean(dataToPlot1st{k,yi}(iPeak,iPlot));
                        tPKyaw = abs(mean(dataToPlot2nd{k,yi}(iPeak,iPlot))-mean(dataToPlot1st{k,yi}(iPeak,iPlot)));
                        trIndex{n,m}=(PKyf-PKyaw)./(tPKyf+tPKyaw);
                end

            case 'Prog-Sideslip'
                switch regions{k}
                    case {'H2axon'}
                        PKss = mean(dataToPlot2nd{k,ssi}(iPeak,iPlot))-mean(dataToPlot1st{k,ssi}(iPeak,iPlot));
                        tPKss = abs(mean(dataToPlot2nd{k,ssi}(iPeak,iPlot))-mean(dataToPlot1st{k,ssi}(iPeak,iPlot)));
                        PKprog = mean(dataToPlot2nd{k,pri}(iPeak,iPlot))-mean(dataToPlot1st{k,pri}(iPeak,iPlot));
                        tPKprog = abs(mean(dataToPlot2nd{k,pri}(iPeak,iPlot))-mean(dataToPlot1st{k,pri}(iPeak,iPlot)));
                        trIndex{n,m}=(PKprog-PKss)./(tPKprog+tPKss);
                    otherwise
                        PKss = mean(dataToPlot1st{k,ssi}(iPeak,iPlot))-mean(dataToPlot2nd{k,ssi}(iPeak,iPlot));
                        tPKss = abs(mean(dataToPlot1st{k,ssi}(iPeak,iPlot))-mean(dataToPlot2nd{k,ssi}(iPeak,iPlot)));
                        PKprog = mean(dataToPlot1st{k,pri}(iPeak,iPlot))-mean(dataToPlot2nd{k,pri}(iPeak,iPlot));
                        tPKprog = abs(mean(dataToPlot1st{k,pri}(iPeak,iPlot))-mean(dataToPlot2nd{k,pri}(iPeak,iPlot)));
                        trIndex{n,m}=(PKprog-PKss)./(tPKprog+tPKss);
                end

            case 'ProgSide45-Sideslip'
                switch regions{k}
                    case {'H2axon'}
                        PKprog45 = mean(dataToPlot2nd{k,ps45i}(iPeak,iPlot))-mean(dataToPlot1st{k,ps45i}(iPeak,iPlot));
                        tPKprog45 = abs(mean(dataToPlot2nd{k,ps45i}(iPeak,iPlot))-mean(dataToPlot1st{k,ps45i}(iPeak,iPlot)));
                        PKss = mean(dataToPlot2nd{k,ssi}(iPeak,iPlot))-mean(dataToPlot1st{k,ssi}(iPeak,iPlot));
                        tPKss = abs(mean(dataToPlot2nd{k,ssi}(iPeak,iPlot))-mean(dataToPlot1st{k,ssi}(iPeak,iPlot)));
                        trIndex{n,m}=(PKprog45-PKss)./(tPKprog45+tPKss);
                    case {'H2rnaxon'}
                        PKprog45 = mean(dataToPlot1st{k,ps45i}(iPeak,iPlot))-mean(dataToPlot2nd{k,ps45i}(iPeak,iPlot));
                        tPKprog45 = abs(mean(dataToPlot1st{k,ps45i}(iPeak,iPlot))-mean(dataToPlot2nd{k,ps45i}(iPeak,iPlot)));
                        PKss = mean(dataToPlot1st{k,ssi}(iPeak,iPlot))-mean(dataToPlot2nd{k,ssi}(iPeak,iPlot));
                        tPKss = abs(mean(dataToPlot1st{k,ssi}(iPeak,iPlot))-mean(dataToPlot2nd{k,ssi}(iPeak,iPlot)));
                        trIndex{n,m}=(PKprog45-PKss)./(tPKprog45+tPKss);
                    otherwise
                        PKprog45 = mean(dataToPlot1st{k,ps45i}(iPeak,iPlot))-mean(dataToPlot2nd{k,ps45i}(iPeak,iPlot));
                        tPKprog45 = abs(mean(dataToPlot1st{k,ps45i}(iPeak,iPlot))-mean(dataToPlot2nd{k,ps45i}(iPeak,iPlot)));
                        PKss = mean(dataToPlot1st{k,ssi}(iPeak,iPlot))-mean(dataToPlot2nd{k,ssi}(iPeak,iPlot));
                        tPKss = abs(mean(dataToPlot1st{k,ssi}(iPeak,iPlot))-mean(dataToPlot2nd{k,ssi}(iPeak,iPlot)));
                        trIndex{n,m}=(PKprog45-PKss)./(tPKprog45+tPKss);
                end
        end
        h(n,m) = notBoxPlot(trIndex{n,m},x(n,m));hold on
        set(h(n,m).data,'MarkerFaceColor',cReg(k,:))
        % save the maximum y values for significance plot alignment
        maxbl = nanmax(trIndex{n,m});
        if isempty(maxbl)
            maxbl = 0;
        end
        my(n) = max(my(n),maxbl);
        trIndex{n,m}(isnan(trIndex{n,m})) = [];
    end
    % Perform Mann-Whitney-Wilcoxon non parametric test and plot
    % significance indicators if p<0.05
    nck = nchoosek(1:m,2);
    %This is all done just to help with the plots
    nck = [nck; nck(end,:)];
    nph = size(nck,1) - 1; % number of comparisons made (used for adjusting the
    % minimum P value required for significance - a.k.a. Bonferroni correction)
    sameComparison = 1;
    for ii = 1 : size(nck,1) - 1
        if nck(ii+1) == nck(ii)
            sameComparison = sameComparison + 1;
        else
            sameComparison = 1;
        end
    end
    for ii = 1 : size(nck,1) - 1
        stats{ii} = mwwtest(trIndex{n,nck(ii,1)},trIndex{n,nck(ii,2)},0);
        % disp([compNames{n},'-',regions{neuronsToPlot(nck(ii,1))},'-',...
        %     regions{neuronsToPlot(nck(ii,2))},'-',num2str(stats{ii}.p(2))])
        if stats{ii}.p < 0.05/nph
            % uncomment below if you want to display the p-values per pair
            disp([compNames{n},'-',regions{neuronsToPlot(nck(ii,1))},'-',...
                regions{neuronsToPlot(nck(ii,2))},'-',num2str(stats{ii}.p(2))])
            sigline([x(n,nck(ii,1)),x(n,nck(ii,2))],[],my(n)+my(n)*0.05*ii);
        end
    end
end
limsX = xlim;
line(limsX,[0,0],'LineStyle','--','Color','k','LineWidth',1)
set(gca,'XTick',x(:,end)','XTickLabel',compNames,'FontSize',18,...
    'YLim',limsY)
set(gca,'FontSize',14,'YLim',limsY)
ylabel('Discrimination Index (peak)')

for m=1:length(neuronsToPlot)
    hold on
    l(m)=plot(nan,nan,'Color',cReg(neuronsToPlot(m),:),'LineWidth',4);
end
set(gca,'XLim',limsX)
lh=legend(l,regions(neuronsToPlot),'Box','off','FontSize',24,'Location',...
    'Best');
clearvars l lh h

%% Figure 5b-e
% --------------------------- User Input ----------------------------------
% Which compound OF stimuli to plot. Defined by their index in the Stimuli
% cell array. Plotted in the same order. Leave blank to plot all stimuli.
% stimToPlot = [];

% stimToPlot = find(contains(Stimuli,{'Trans'})&contains(Stimuli,{'5px'}));
% stimToPlot(contains(Stimuli(stimToPlot),{'Trans1Rot0'})) = [];
% stimToPlot(contains(Stimuli(stimToPlot),{'Trans4-28'})) = [];

stimToPlot=find(contains(Stimuli,{'Trans0Rot1'})&contains(Stimuli,{'5px'}));

% which neurons to plot. Defined by their index in the regions cell array.
% Plotted in the same order. Leave blank to plot all neurons.

neuronsToPlot = [6,7,9,10,1,2]; % Figure 5b neurons

% Which calcium indicator(s) to plot. Defined by the genotype field in OF
% struct array. Leave blank to plot everything. Note that it is possible to
% pool data from multiple indicators (e.g. gcampsToPlot = {GC6f','GC7f'})
% Accepted arguments: {'GC6f','GC7f','sytGC7f'} Default: gcampsToPlot = {};
gcampsToPlot = {'GC6f'};

% Which prep types to plot. Defined by the prepType field in OF struct
% array. Leave blank to plot all experiments.
% Accepted arguments: ['walking','immobilized'] Default: prepToPlot = '';
prepToPlot = '';

% Neurons from which hemisphere to plot. Defined by the somaSide field in
% OF struct array. Type 'Both' to plot all neurons. In that case, responses
% from Right hemisphere will be flipped to match the Left hemisphere neural
% responses and plotted together.
% Accepted arguments: ['L','R','Both'] Default: sideToPlot = 'Both';
sideToPlot = 'Both';

% specify normalization type. Type 'Raw' for no normalization.
% 'Max' will normalize the data per ROI to the maximum response recorded in
% that ROI (see variable maxes).
% 'Zscore' will zscore the data per ROI based on the mean and standard
% deviation of the entire recording session done on that ROI (see variables
% mu and sigma). 'Znorm' normalizes Z-score data by the max Z-score
% observed per fly per neuron
% Accepted arguments: ['Max','Zscore','Raw','Znorm'] Default: 'Raw';
normType = 'Zscore';

% Specify plot type. Type 'Line' for plotting each fly with a different
% colored line together with the mean trace. 'Shaded' will plot the std as
% a shaded area around the mean trace instead.
% Accepted arguments: ['Line','Shaded']
plotType = 'Shaded';

switch normType
    case 'Max'
        limsY = [-0.5 1.1]; % Y-axis limits
    case 'Raw'
        limsY = [-3 10]; % Y-axis limits
    case 'Zscore'
        switch plotType
            case 'Line'
                limsY = [-5 30]; % Y-axis limits
            case 'Shaded'
                limsY = [-2.5 18]; % Y-axis limits
        end
    case 'Znorm'
        limsY = [-0.25 1]; % Y-axis limits
end
%----------------------------end of user input-----------------------------

% Define the subset of data that will be plotted. This will be used in the
% plotSelectDFF function to plot only the specified subset of data as well
% as to calculate and display number of flies & ROIs being plotted
if isempty(stimToPlot)
    stimToPlot = 1:length(Stimuli);
end
if isempty(neuronsToPlot)
    neuronsToPlot = 1:length(regions);
end

prepIndices = cell(length(regions),1);
gcIndices = cell(length(regions),1);
sideIndices = cell(length(regions),1);
prepIndicesFly = cell(length(regions),1);
gcIndicesFly = cell(length(regions),1);
sideIndicesFly = cell(length(regions),1);

if ~isempty(prepToPlot)
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        if isempty(prepIndices{r})
            prepIndices{r} = ...
                contains(PrepTypePerROI{r},prepToPlot);
            prepIndicesFly{r} = ...
                contains(PrepType{r},prepToPlot);
        else
            prepIndices{r} = (prepIndices{r} | ...
                contains(PrepTypePerROI{r},prepToPlot));
            prepIndicesFly{r} = (prepIndicesFly{r} | ...
                contains(PrepType{r},prepToPlot));
        end
    end
else
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        prepIndices{r} = true(size(PrepTypePerROI{r}));
        prepIndicesFly{r} = true(size(PrepType{r}));
    end
end
if ~isempty(gcampsToPlot)
    for p = 1:length(gcampsToPlot) %for each selected GCaMP type
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            if isempty(gcIndices{r})
                gcIndices{r} = ...
                    ismember(GCaMPTypePerROI{r},gcampsToPlot{p});
                gcIndicesFly{r} = ...
                    ismember(GCaMPType{r},gcampsToPlot{p});
            else
                gcIndices{r} = (gcIndices{r} | ...
                    ismember(GCaMPTypePerROI{r},gcampsToPlot{p}));
                gcIndicesFly{r} = (gcIndicesFly{r} | ...
                    ismember(GCaMPType{r},gcampsToPlot{p}));
            end
        end
    end
else
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        gcIndices{r} = true(size(GCaMPTypePerROI{r}));
        gcIndicesFly{r} = true(size(GCaMPType{r}));
    end
end

switch sideToPlot
    case 'Both'
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            sideIndices{r} = true(size(SomaSidePerROI{r}));
            sideIndicesFly{r} = true(size(SomaSide{r}));
        end
        switch normType
            case 'Zscore'
                dataToPlot = OFRes_1st_Z_R;
                aucToPlot = Auc_1st_bs_Z_R;
            case 'Max'
                dataToPlot = STResP_normM;
            case 'Raw'
                dataToPlot = OFRes_1st;
            case 'Znorm'
                dataToPlot = OFRes_1st_Znorm;
            otherwise
                error(['Normalization type not specified correctly. ',...
                    'Use normType=[''Max'',''Zscore'',''Znorm'', ',...
                    'or ''Raw'']'])
        end
    case {'L','R'}
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            if isempty(sideIndices{r})
                sideIndices{r} = ...
                    contains(SomaSidePerROI{r},sideToPlot);
                sideIndicesFly{r} = ...
                    contains(SomaSide{r},sideToPlot);
            else
                sideIndices{r} = (sideIndices{r} | ...
                    contains(SomaSidePerROI{r},sideToPlot));
                sideIndicesFly{r} = (sideIndicesFly{r} | ...
                    contains(SomaSide{r},sideToPlot));
            end
        end
        switch normType
            case 'Zscore'
                dataToPlot = OFRes_Raw_1st_Z;
            case 'Max'
                dataToPlot = OFRes_Raw_1st_normM;
            case 'Raw'
                dataToPlot1st = OFRes_Raw_1st;
            otherwise
                error(['Normalization type not specified correctly. ',...
                    'Use normType=[''Max'',''Zscore'', ',...
                    'or ''Raw'']'])
        end
    otherwise
        error(['Plotting side not specified correctly. ',...
            'Use sideToPlot = [''L'',''R'' or ''Both'']'])
end

f=plotSelectDFFCompound(dataToPlot,neuronsToPlot,...
    stimToPlot,minFrameRate,cReg,prepIndices,prepIndicesFly,...
    gcIndices,gcIndicesFly,sideIndices,sideIndicesFly,Stimuli,...
    regions,limsY,plotType);
% -------------------------------------------------------------------------
% Figure 5c
stimToPlot=find(contains(Stimuli,{'Trans0-71','Trans1-42','Trans2-14'})&...
    contains(Stimuli,{'5px'}));
stimToPlot(contains(Stimuli(stimToPlot),{'Rot0_CW'})) = [];
neuronsToPlot = [6,7];

f=plotSelectDFFCompound(dataToPlot,neuronsToPlot,...
    stimToPlot,minFrameRate,cReg,prepIndices,prepIndicesFly,...
    gcIndices,gcIndicesFly,sideIndices,sideIndicesFly,Stimuli,...
    regions,limsY,plotType);
% -------------------------------------------------------------------------
% Figure 5d
neuronsToPlot = [9,10];
f=plotSelectDFFCompound(dataToPlot,neuronsToPlot,...
    stimToPlot,minFrameRate,cReg,prepIndices,prepIndicesFly,...
    gcIndices,gcIndicesFly,sideIndices,sideIndicesFly,Stimuli,...
    regions,limsY,plotType);
% -------------------------------------------------------------------------
% Figure 5e
stimToPlot=find(contains(Stimuli,{'Trans0-71','Trans1-42','Trans2-14'})&...
    contains(Stimuli,{'5px'}));
stimToPlot(contains(Stimuli(stimToPlot),{'Rot0_CW','Rot1_CCW'})) = [];
neuronsToPlot = [1,2,6];
f=plotSelectDFFCompound(dataToPlot,neuronsToPlot,...
    stimToPlot,minFrameRate,cReg,prepIndices,prepIndicesFly,...
    gcIndices,gcIndicesFly,sideIndices,sideIndicesFly,Stimuli,...
    regions,limsY,plotType);
% -------------------------------------------------------------------------
%% Figure 5f

stimToPlot = find(contains(Stimuli,{'Trans'})&contains(Stimuli,{'5px'}));
stimToPlot(contains(Stimuli(stimToPlot),{'Trans1Rot0','Trans4-28'})) = [];
stimToPlot(contains(Stimuli(stimToPlot),{'Rot1','Flicker','_CCW'})) = [];

% which neurons to plot. Defined by their index in the regions cell array.
% Plotted in the same order. Leave blank to plot all neurons.

neuronsToPlot = [6,7,9,10,1,2];

% Which calcium indicator(s) to plot. Defined by the genotype field in OF
% struct array. Leave blank to plot everything. Note that it is possible to
% pool data from multiple indicators (e.g. gcampsToPlot = {GC6f','GC7f'})
% Accepted arguments: {'GC6f','GC7f','sytGC7f'} Default: gcampsToPlot = {};
gcampsToPlot = {'GC6f'};

% Which prep types to plot. Defined by the prepType field in OF struct
% array. Leave blank to plot all experiments.
% Accepted arguments: ['walking','immobilized'] Default: prepToPlot = '';
prepToPlot = '';

% Neurons from which hemisphere to plot. Defined by the somaSide field in
% OF struct array. Type 'Both' to plot all neurons. In that case, responses
% from Right hemisphere will be flipped to match the Left hemisphere neural
% responses and plotted together.
% Accepted arguments: ['L','R','Both'] Default: sideToPlot = 'Both';
sideToPlot = 'Both';

% specify normalization type. Type 'Raw' for no normalization.
% 'Max' will normalize the data per ROI to the maximum response recorded in
% that ROI (see variable maxes).
% 'Zscore' will zscore the data per ROI based on the mean and standard
% deviation of the entire recording session done on that ROI (see variables
% mu and sigma). 'Znorm' normalizes Z-score data by the max Z-score
% observed per fly per neuron
% Accepted arguments: ['Max','Zscore','Raw','Znorm'] Default: 'Raw';
normType = 'Zscore';

%----------------------------end of user input-----------------------------

peakbin = 500; % time bin for the averaged peak response (ms)
H2Flip = 0;
iH2 = find(contains(regions,'H2axon')); % index of H2
inotH2 = find(~contains(regions,'H2axon'));
% index for recordings that will come from the right hemisphere
[~,loc] = ismember(StimuliFlipped,Stimuli);

clearvars l lh h stats

% Define the subset of data that will be plotted. This will be used in the
% plotSelectAUC function to plot only the specified subset of data as well
% as to calculate and display number of flies & ROIs being plotted
if isempty(stimToPlot)
    stimToPlot = 1:length(Stimuli);
end
if isempty(neuronsToPlot)
    neuronsToPlot = 1:length(regions);
end

prepIndices = cell(length(regions),1);
gcIndices = cell(length(regions),1);
sideIndices = cell(length(regions),1);
prepIndicesFly = cell(length(regions),1);
gcIndicesFly = cell(length(regions),1);
sideIndicesFly = cell(length(regions),1);

if ~isempty(prepToPlot)
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        if isempty(prepIndices{r})
            prepIndices{r} = ...
                contains(PrepTypePerROI{r},prepToPlot);
            prepIndicesFly{r} = ...
                contains(PrepType{r},prepToPlot);
        else
            prepIndices{r} = (prepIndices{r} | ...
                contains(PrepTypePerROI{r},prepToPlot));
            prepIndicesFly{r} = (prepIndicesFly{r} | ...
                contains(PrepType{r},prepToPlot));
        end
    end
else
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        prepIndices{r} = true(size(PrepTypePerROI{r}));
        prepIndicesFly{r} = true(size(PrepType{r}));
    end
end
if ~isempty(gcampsToPlot)
    for p = 1:length(gcampsToPlot) %for each selected GCaMP type
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            if isempty(gcIndices{r})
                gcIndices{r} = ...
                    ismember(GCaMPTypePerROI{r},gcampsToPlot{p});
                gcIndicesFly{r} = ...
                    ismember(GCaMPType{r},gcampsToPlot{p});
            else
                gcIndices{r} = (gcIndices{r} | ...
                    ismember(GCaMPTypePerROI{r},gcampsToPlot{p}));
                gcIndicesFly{r} = (gcIndicesFly{r} | ...
                    ismember(GCaMPType{r},gcampsToPlot{p}));
            end
        end
    end
else
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        gcIndices{r} = true(size(GCaMPTypePerROI{r}));
        gcIndicesFly{r} = true(size(GCaMPType{r}));
    end
end
switch sideToPlot
    case 'Both'
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            sideIndices{r} = true(size(SomaSidePerROI{r}));
            sideIndicesFly{r} = true(size(SomaSide{r}));
        end
        switch normType
            case 'Zscore'
                aucToPlot = Auc_1st_bs_Z_R;
                dataToPlot = OFRes_1st_Z_R;
            case 'Max'
                dataToPlot = Auc_1st_normM;
            case 'Raw'
                dataToPlot1st = Auc_1st_bs;
                dataToPlot2nd = Auc_2nd_bs;
            otherwise
                error(['Normalization type not specified correctly. ',...
                    'Use normType=[''Max'',''Zscore'', ',...
                    'or ''Raw'']'])
        end
    case {'L','R'}
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            if isempty(sideIndices{r})
                sideIndices{r} = ...
                    contains(SomaSidePerROI{r},sideToPlot);
                sideIndicesFly{r} = ...
                    contains(SomaSide{r},sideToPlot);
            else
                sideIndices{r} = (sideIndices{r} | ...
                    contains(SomaSidePerROI{r},sideToPlot));
                sideIndicesFly{r} = (sideIndicesFly{r} | ...
                    contains(SomaSide{r},sideToPlot));
            end
        end
        switch normType
            case 'Zscore'
                dataToPlot1st = Auc_1st_bs_Z;
                dataToPlot2nd = Auc_2nd_bs_Z;
            case 'Max'
                dataToPlot1st = Auc_1st_normM;
                dataToPlot2nd = Auc_2nd_normM;
            case 'Raw'
                dataToPlot1st = Auc_1st_bs;
                dataToPlot2nd = Auc_2nd_bs;
            otherwise
                error(['Normalization type not specified correctly. ',...
                    'Use normType=[''Max'',''Zscore'', ',...
                    'or ''Raw'']'])
        end
    otherwise
        error(['Plotting side not specified correctly. ',...
            'Use sideToPlot = [''L'',''R'' or ''Both'']'])
end

frameNum = size(dataToPlot{1,1},1);
TimeAxisPD = linspace(0,frameNum/(minFrameRate/1000),frameNum);
limsY = [-2.5 20];
iPeak = find(TimeAxisPD>4-peakbin/1000,1,'first'):...
    find(TimeAxisPD<4,1,'last');

sp=0.5;
f=figure('WindowState','maximized');
StimXLabel = {'Low','Intermediate','Fast'};
x = zeros(size(neuronsToPlot,2),size(stimToPlotCW,2));
DSI = cell(size(neuronsToPlot,2),size(stimToPlotCW,2));
DSR = cell(size(neuronsToPlot,2),size(stimToPlotCW,2));
NPD = cell(size(neuronsToPlot,2),size(stimToPlotCW,2));
my = zeros(1,size(neuronsToPlot,2));

for n = 1:length(neuronsToPlot)
    s = neuronsToPlot(n);
    iPlot = prepIndices{s} & gcIndices{s} & sideIndices{s};
    for k = 1:length(stimToPlotCW)
        mCW = stimToPlotCW(k);
        mCCW = stimToPlotCCW(k);
        switch regions{s}
            case {'H2axon'}
                rPD = mean(dataToPlot{s,mCCW}(iPeak,iPlot));
                rNPD = mean(dataToPlot{s,mCW}(iPeak,iPlot));
            otherwise
                rPD = mean(dataToPlot{s,mCW}(iPeak,iPlot));
                rNPD = mean(dataToPlot{s,mCCW}(iPeak,iPlot));
        end
        DSI{n,k}=(rPD-rNPD)./(abs(rPD)+abs(rNPD));
        DSR{n,k}=(rPD-rNPD);
        NPD{n,k}=(rNPD);
        maxbl = nanmax(DSR{n,k});
        if isempty(maxbl)
            maxbl = 0;
        end
        my(n) = max(my(n),maxbl);
        NPD{n,k}(isnan(NPD{n,k})) = [];
    end
    for k = 1:length(stimToPlot)
        if isempty(stimToPlotCW==stimToPlot(k))
            sti = find(stimToPlotCCW==stimToPlot(k));
        else
            sti = find(stimToPlotCW==stimToPlot(k));
        end
        x(n,k) = (n-1)*length(stimToPlot) + 1 + (k-1)*sp;
        h(n,k) = notBoxPlot(NPD{n,sti},x(n,k));hold on
        set(h(n,k).data,'MarkerFaceColor',cReg(s,:))
    end
end
limsX = xlim;
line(limsX,[0,0],'LineStyle','--','Color','k','LineWidth',1)
set(gca,'XTick',unique(x),'XTickLabel',StimXLabel,...
    'FontSize',14,'TickLabelInterpreter','none',...
    'XTickLabelRotation',45)
set(gca,'YLim',limsY)
ylabel('Straight translation Response (z-score)')

for m=1:length(neuronsToPlot)
    hold on
    l(m)=plot(nan,nan,'Color',cReg(neuronsToPlot(m),:),'LineWidth',4);
end
set(gca,'XLim',limsX)
lh=legend(l,regions(neuronsToPlot),'Box','off','FontSize',20,...
    'Location','Best','Interpreter','none');

%% Figure 5g

% Which compound OF stimuli to plot. Defined by their index in the Stimuli
% cell array. Plotted in the same order. Leave blank to plot all stimuli.

stimToPlot = find(contains(Stimuli,{'Trans1-42Rot1'})&contains(Stimuli,{'5px'}));

% which neurons to plot. Defined by their index in the regions cell array.
% Plotted in the same order. Leave blank to plot all neurons.

neuronsToPlot = [6,7,9,10,1,2];

% Which calcium indicator(s) to plot. Defined by the genotype field in OF
% struct array. Leave blank to plot everything. Note that it is possible to
% pool data from multiple indicators (e.g. gcampsToPlot = {GC6f','GC7f'})
% Accepted arguments: {'GC6f','GC7f','sytGC7f'} Default: gcampsToPlot = {};
gcampsToPlot = {'GC6f'};

% Which prep types to plot. Defined by the prepType field in OF struct
% array. Leave blank to plot all experiments.
% Accepted arguments: ['walking','immobilized'] Default: prepToPlot = '';
prepToPlot = '';

% Neurons from which hemisphere to plot. Defined by the somaSide field in
% OF struct array. Type 'Both' to plot all neurons. In that case, responses
% from Right hemisphere will be flipped to match the Left hemisphere neural
% responses and plotted together.
% Accepted arguments: ['L','R','Both'] Default: sideToPlot = 'Both';
sideToPlot = 'Both';

% specify normalization type. Type 'Raw' for no normalization.
% 'Max' will normalize the data per ROI to the maximum response recorded in
% that ROI (see variable maxes).
% 'Zscore' will zscore the data per ROI based on the mean and standard
% deviation of the entire recording session done on that ROI (see variables
% mu and sigma). 'Znorm' normalizes Z-score data by the max Z-score
% observed per fly per neuron
% Accepted arguments: ['Max','Zscore','Raw','Znorm'] Default: 'Raw';
normType = 'Zscore';

%----------------------------end of user input-----------------------------

peakbin = 500; % time bin for the averaged peak response (ms)
H2Flip = 0;
iH2 = find(contains(regions,'H2axon')); % index of H2
inotH2 = find(~contains(regions,'H2axon'));
% index for recordings that will come from the right hemisphere
[~,loc] = ismember(StimuliFlipped,Stimuli);

clearvars l lh h stats

% Define the subset of data that will be plotted. This will be used in the
% plotSelectAUC function to plot only the specified subset of data as well
% as to calculate and display number of flies & ROIs being plotted
if isempty(stimToPlot)
    stimToPlot = 1:length(Stimuli);
end
if isempty(neuronsToPlot)
    neuronsToPlot = 1:length(regions);
end

prepIndices = cell(length(regions),1);
gcIndices = cell(length(regions),1);
sideIndices = cell(length(regions),1);
prepIndicesFly = cell(length(regions),1);
gcIndicesFly = cell(length(regions),1);
sideIndicesFly = cell(length(regions),1);

if ~isempty(prepToPlot)
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        if isempty(prepIndices{r})
            prepIndices{r} = ...
                contains(PrepTypePerROI{r},prepToPlot);
            prepIndicesFly{r} = ...
                contains(PrepType{r},prepToPlot);
        else
            prepIndices{r} = (prepIndices{r} | ...
                contains(PrepTypePerROI{r},prepToPlot));
            prepIndicesFly{r} = (prepIndicesFly{r} | ...
                contains(PrepType{r},prepToPlot));
        end
    end
else
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        prepIndices{r} = true(size(PrepTypePerROI{r}));
        prepIndicesFly{r} = true(size(PrepType{r}));
    end
end
if ~isempty(gcampsToPlot)
    for p = 1:length(gcampsToPlot) %for each selected GCaMP type
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            if isempty(gcIndices{r})
                gcIndices{r} = ...
                    ismember(GCaMPTypePerROI{r},gcampsToPlot{p});
                gcIndicesFly{r} = ...
                    ismember(GCaMPType{r},gcampsToPlot{p});
            else
                gcIndices{r} = (gcIndices{r} | ...
                    ismember(GCaMPTypePerROI{r},gcampsToPlot{p}));
                gcIndicesFly{r} = (gcIndicesFly{r} | ...
                    ismember(GCaMPType{r},gcampsToPlot{p}));
            end
        end
    end
else
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        gcIndices{r} = true(size(GCaMPTypePerROI{r}));
        gcIndicesFly{r} = true(size(GCaMPType{r}));
    end
end
switch sideToPlot
    case 'Both'
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            sideIndices{r} = true(size(SomaSidePerROI{r}));
            sideIndicesFly{r} = true(size(SomaSide{r}));
        end
        switch normType
            case 'Zscore'
                aucToPlot = Auc_1st_bs_Z_R;
                dataToPlot = OFRes_1st_Z_R;
            case 'Max'
                dataToPlot = Auc_1st_normM;
            case 'Raw'
                dataToPlot1st = Auc_1st_bs;
                dataToPlot2nd = Auc_2nd_bs;
            otherwise
                error(['Normalization type not specified correctly. ',...
                    'Use normType=[''Max'',''Zscore'', ',...
                    'or ''Raw'']'])
        end
    case {'L','R'}
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            if isempty(sideIndices{r})
                sideIndices{r} = ...
                    contains(SomaSidePerROI{r},sideToPlot);
                sideIndicesFly{r} = ...
                    contains(SomaSide{r},sideToPlot);
            else
                sideIndices{r} = (sideIndices{r} | ...
                    contains(SomaSidePerROI{r},sideToPlot));
                sideIndicesFly{r} = (sideIndicesFly{r} | ...
                    contains(SomaSide{r},sideToPlot));
            end
        end
        switch normType
            case 'Zscore'
                dataToPlot1st = Auc_1st_bs_Z;
                dataToPlot2nd = Auc_2nd_bs_Z;
            case 'Max'
                dataToPlot1st = Auc_1st_normM;
                dataToPlot2nd = Auc_2nd_normM;
            case 'Raw'
                dataToPlot1st = Auc_1st_bs;
                dataToPlot2nd = Auc_2nd_bs;
            otherwise
                error(['Normalization type not specified correctly. ',...
                    'Use normType=[''Max'',''Zscore'', ',...
                    'or ''Raw'']'])
        end
    otherwise
        error(['Plotting side not specified correctly. ',...
            'Use sideToPlot = [''L'',''R'' or ''Both'']'])
end

frameNum = size(dataToPlot{1,1},1);
TimeAxisPD = linspace(0,frameNum/(minFrameRate/1000),frameNum);
limsY = [-2.5 25];
iPeak = find(TimeAxisPD>4-peakbin/1000,1,'first'):...
    find(TimeAxisPD<4,1,'last');

sp=0.5;
f=figure('WindowState','maximized');
StimXLabel = {'Left','Right'};
x = zeros(size(neuronsToPlot,2),size(stimToPlotCW,2));
DSI = cell(size(neuronsToPlot,2),size(stimToPlotCW,2));
DSR = cell(size(neuronsToPlot,2),size(stimToPlotCW,2));
DSRnormM = cell(size(neuronsToPlot,2),size(stimToPlotCW,1));
DSRnormMA = cell(size(neuronsToPlot,2),size(stimToPlotCW,1));
NPD = cell(size(neuronsToPlot,2),size(stimToPlotCW,2));
PD = cell(size(neuronsToPlot,2),size(stimToPlotCW,2));
my = zeros(1,size(neuronsToPlot,2));

for n = 1:length(neuronsToPlot)
    s = neuronsToPlot(n);
    iPlot = prepIndices{s} & gcIndices{s} & sideIndices{s};
    for k = 1:length(stimToPlotCW)
        mCW = stimToPlotCW(k);
        mCCW = stimToPlotCCW(k);
        switch regions{s}
            case {'H2axon'}
                rPD = mean(dataToPlot{s,mCCW}(iPeak,iPlot));
                rNPD = mean(dataToPlot{s,mCW}(iPeak,iPlot));
            otherwise
                rPD = mean(dataToPlot{s,mCW}(iPeak,iPlot));
                rNPD = mean(dataToPlot{s,mCCW}(iPeak,iPlot));
        end
        DSI{n,k}=(rPD-rNPD)./(abs(rPD)+abs(rNPD));
        DSR{n,k}=(rPD-rNPD);
        PD{n,k}=(rPD);
        NPD{n,k}=(rNPD);
        maxbl = nanmax(DSR{n,k});
        if isempty(maxbl)
            maxbl = 0;
        end
        my(n) = max(my(n),maxbl);
        NPD{n,k}(isnan(NPD{n,k})) = [];
        PD{n,k}(isnan(PD{n,k})) = [];
    end
    for k = 1:length(stimToPlot)
        if sum(stimToPlotCW==stimToPlot(k))
            sti = find(stimToPlotCW==stimToPlot(k));
            x(n,k) = (n-1)*length(stimToPlot) + 1 + (k-1)*sp;
            h(n,k) = notBoxPlot(NPD{n,sti},x(n,k));hold on
            set(h(n,k).data,'MarkerFaceColor',cReg(s,:))
        else
            sti = find(stimToPlotCCW==stimToPlot(k));
            x(n,k) = (n-1)*length(stimToPlot) + 1 + (k-1)*sp;
            h(n,k) = notBoxPlot(PD{n,sti},x(n,k));hold on
            set(h(n,k).data,'MarkerFaceColor',cReg(s,:))
        end
    end
end
limsX = xlim;
line(limsX,[0,0],'LineStyle','--','Color','k','LineWidth',1)
set(gca,'XTick',unique(x),'XTickLabel',StimXLabel,...
    'FontSize',14,'TickLabelInterpreter','none',...
    'XTickLabelRotation',45)
set(gca,'YLim',limsY)
ylabel('Deviated translation Response (z-score)')

for m=1:length(neuronsToPlot)
    hold on
    l(m)=plot(nan,nan,'Color',cReg(neuronsToPlot(m),:),'LineWidth',4);
end
set(gca,'XLim',limsX)
lh=legend(l,regions(neuronsToPlot),'Box','off','FontSize',20,...
    'Location','Best','Interpreter','none');

%% Figure 5h

% Which compound OF stimuli to plot. Defined by their index in the Stimuli
% cell array. Plotted in the same order. Leave blank to plot all stimuli.

stimToPlot = find(contains(Stimuli,{'Trans1-42Rot1'})&contains(Stimuli,{'5px'}));

% which neurons to plot. Defined by their index in the regions cell array.
% Plotted in the same order. Leave blank to plot all neurons.

neuronsToPlot = [6,7,9,10,1,2];

% Which calcium indicator(s) to plot. Defined by the genotype field in OF
% struct array. Leave blank to plot everything. Note that it is possible to
% pool data from multiple indicators (e.g. gcampsToPlot = {GC6f','GC7f'})
% Accepted arguments: {'GC6f','GC7f','sytGC7f'} Default: gcampsToPlot = {};
gcampsToPlot = {'GC6f'};

% Which prep types to plot. Defined by the prepType field in OF struct
% array. Leave blank to plot all experiments.
% Accepted arguments: ['walking','immobilized'] Default: prepToPlot = '';
prepToPlot = '';

% Neurons from which hemisphere to plot. Defined by the somaSide field in
% OF struct array. Type 'Both' to plot all neurons. In that case, responses
% from Right hemisphere will be flipped to match the Left hemisphere neural
% responses and plotted together.
% Accepted arguments: ['L','R','Both'] Default: sideToPlot = 'Both';
sideToPlot = 'Both';

% specify normalization type. Type 'Raw' for no normalization.
% 'Max' will normalize the data per ROI to the maximum response recorded in
% that ROI (see variable maxes).
% 'Zscore' will zscore the data per ROI based on the mean and standard
% deviation of the entire recording session done on that ROI (see variables
% mu and sigma). 'Znorm' normalizes Z-score data by the max Z-score
% observed per fly per neuron
% Accepted arguments: ['Max','Zscore','Raw','Znorm'] Default: 'Raw';
normType = 'Zscore';

%----------------------------end of user input-----------------------------

peakbin = 500; % time bin for the averaged peak response (ms)

clearvars l lh h stats

% Define the subset of data that will be plotted. This will be used in the
% plotSelectAUC function to plot only the specified subset of data as well
% as to calculate and display number of flies & ROIs being plotted
if isempty(stimToPlot)
    stimToPlot = 1:length(Stimuli);
end
if isempty(neuronsToPlot)
    neuronsToPlot = 1:length(regions);
end

prepIndices = cell(length(regions),1);
gcIndices = cell(length(regions),1);
sideIndices = cell(length(regions),1);
prepIndicesFly = cell(length(regions),1);
gcIndicesFly = cell(length(regions),1);
sideIndicesFly = cell(length(regions),1);

if ~isempty(prepToPlot)
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        if isempty(prepIndices{r})
            prepIndices{r} = ...
                contains(PrepTypePerROI{r},prepToPlot);
            prepIndicesFly{r} = ...
                contains(PrepType{r},prepToPlot);
        else
            prepIndices{r} = (prepIndices{r} | ...
                contains(PrepTypePerROI{r},prepToPlot));
            prepIndicesFly{r} = (prepIndicesFly{r} | ...
                contains(PrepType{r},prepToPlot));
        end
    end
else
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        prepIndices{r} = true(size(PrepTypePerROI{r}));
        prepIndicesFly{r} = true(size(PrepType{r}));
    end
end
if ~isempty(gcampsToPlot)
    for p = 1:length(gcampsToPlot) %for each selected GCaMP type
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            if isempty(gcIndices{r})
                gcIndices{r} = ...
                    ismember(GCaMPTypePerROI{r},gcampsToPlot{p});
                gcIndicesFly{r} = ...
                    ismember(GCaMPType{r},gcampsToPlot{p});
            else
                gcIndices{r} = (gcIndices{r} | ...
                    ismember(GCaMPTypePerROI{r},gcampsToPlot{p}));
                gcIndicesFly{r} = (gcIndicesFly{r} | ...
                    ismember(GCaMPType{r},gcampsToPlot{p}));
            end
        end
    end
else
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        gcIndices{r} = true(size(GCaMPTypePerROI{r}));
        gcIndicesFly{r} = true(size(GCaMPType{r}));
    end
end
switch sideToPlot
    case 'Both'
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            sideIndices{r} = true(size(SomaSidePerROI{r}));
            sideIndicesFly{r} = true(size(SomaSide{r}));
        end
        switch normType
            case 'Zscore'
                aucToPlot = Auc_1st_bs_Z_R;
                dataToPlot = OFRes_1st_Z_R;
            case 'Max'
                dataToPlot = Auc_1st_normM;
            case 'Raw'
                dataToPlot1st = Auc_1st_bs;
                dataToPlot2nd = Auc_2nd_bs;
            otherwise
                error(['Normalization type not specified correctly. ',...
                    'Use normType=[''Max'',''Zscore'', ',...
                    'or ''Raw'']'])
        end
    case {'L','R'}
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            if isempty(sideIndices{r})
                sideIndices{r} = ...
                    contains(SomaSidePerROI{r},sideToPlot);
                sideIndicesFly{r} = ...
                    contains(SomaSide{r},sideToPlot);
            else
                sideIndices{r} = (sideIndices{r} | ...
                    contains(SomaSidePerROI{r},sideToPlot));
                sideIndicesFly{r} = (sideIndicesFly{r} | ...
                    contains(SomaSide{r},sideToPlot));
            end
        end
        switch normType
            case 'Zscore'
                dataToPlot1st = Auc_1st_bs_Z;
                dataToPlot2nd = Auc_2nd_bs_Z;
            case 'Max'
                dataToPlot1st = Auc_1st_normM;
                dataToPlot2nd = Auc_2nd_normM;
            case 'Raw'
                dataToPlot1st = Auc_1st_bs;
                dataToPlot2nd = Auc_2nd_bs;
            otherwise
                error(['Normalization type not specified correctly. ',...
                    'Use normType=[''Max'',''Zscore'', ',...
                    'or ''Raw'']'])
        end
    otherwise
        error(['Plotting side not specified correctly. ',...
            'Use sideToPlot = [''L'',''R'' or ''Both'']'])
end

frameNum = size(dataToPlot{1,1},1);
TimeAxisPD = linspace(0,frameNum/(minFrameRate/1000),frameNum);
limsY = [-5 30];
iPeak = find(TimeAxisPD>4-peakbin/1000,1,'first'):...
    find(TimeAxisPD<4,1,'last');
sp=0.5;
f=figure('WindowState','maximized');
x = zeros(size(neuronsToPlot,2),size(stimToPlotCW,2));
DSI = cell(size(neuronsToPlot,2),size(stimToPlotCW,2));
DSR = cell(size(neuronsToPlot,2),size(stimToPlotCW,2));
DSRnormM = cell(size(neuronsToPlot,2),size(stimToPlotCW,1));
DSRnormMA = cell(size(neuronsToPlot,2),size(stimToPlotCW,1));
my = zeros(1,size(neuronsToPlot,2));
StimXLabel = {'zero','low','int.','fast'};
for n = 1:length(neuronsToPlot)
    s = neuronsToPlot(n);
    iPlot = prepIndices{s} & gcIndices{s} & sideIndices{s};
    for k = 1:length(stimToPlotCW)
        mCW = stimToPlotCW(k);
        mCCW = stimToPlotCCW(k);
        switch regions{s}
            case {'H2axon'}
                rPD = mean(dataToPlot{s,mCCW}(iPeak,iPlot));
                rNPD = mean(dataToPlot{s,mCW}(iPeak,iPlot));
            otherwise
                rPD = mean(dataToPlot{s,mCW}(iPeak,iPlot));
                rNPD = mean(dataToPlot{s,mCCW}(iPeak,iPlot));
        end
        DSI{n,k}=(rPD-rNPD)./(abs(rPD)+abs(rNPD));
        DSR{n,k}=(rPD-rNPD);
        maxbl = nanmax(DSR{n,k});
        if isempty(maxbl)
            maxbl = 0;
        end
        my(n) = max(my(n),maxbl);
    end
    maxnorm = zeros(1,size(DSR{n,k},2));
    for i = 1:size(DSR{n,k},2)
        for k = 1:length(stimToPlotCW)
            if i <= size(DSR{n,k},2)
                maxnorm(1,i) = max([maxnorm(1,i), abs(DSR{n,k}(i))]);
            end
        end
    end
    for k = 1:length(stimToPlotCW)
        DSRnormM{n,k}=DSR{n,k}./maxnorm(1,1:size(DSR{n,k},2));
        DSRnormMA{n,k}=abs(DSR{n,k})./maxnorm(1,1:size(DSR{n,k},2));
        DSI{n,k}(isnan(DSI{n,k})) = [];
        DSR{n,k}(isnan(DSR{n,k})) = [];
        DSRnormM{n,k}(isnan(DSRnormM{n,k})) = [];
        DSRnormMA{n,k}(isnan(DSRnormMA{n,k})) = [];
    end
    for k = 1:4
        sti = k+5;
        x(n,k) = (n-1)*length(stimToPlot) + 1 + (k-1)*sp;
        h(n,k) = notBoxPlot(DSR{n,sti},x(n,k));hold on
        set(h(n,k).data,'MarkerFaceColor',cReg(s,:))
    end
end
limsX = xlim;
line(limsX,[0,0],'LineStyle','--','Color','k','LineWidth',1)
set(gca,'XTick',unique(x),'XTickLabel',StimXLabel,...
    'FontSize',14,'TickLabelInterpreter','none',...
    'XTickLabelRotation',45)
set(gca,'YLim',limsY)
ylabel('Direction Selective Response (z-score)')

for m=1:length(neuronsToPlot)
    hold on
    l(m)=plot(nan,nan,'Color',cReg(neuronsToPlot(m),:),'LineWidth',4);
end
set(gca,'XLim',limsX)
lh=legend(l,regions(neuronsToPlot),'Box','off','FontSize',20,...
    'Location','Best','Interpreter','none');

%% Figure 6g

% which neurons to plot. Defined by their index in the regions cell array.
% Plotted in the same order. Leave blank to plot all neurons.

neuronsToPlot = [9,10];

% Which calcium indicator(s) to plot. Defined by the genotype field in OF
% struct array. Leave blank to plot everything. Note that it is possible to
% pool data from multiple indicators (e.g. gcampsToPlot = {GC6f','GC7f'})
% Accepted arguments: {'GC6f','GC7f','sytGC7f'} Default: gcampsToPlot = {};
gcampsToPlot = {'GC6f'};

% Which prep types to plot. Defined by the prepType field in OF struct
% array. Leave blank to plot all experiments.
% Accepted arguments: ['walking','immobilized'] Default: prepToPlot = '';
prepToPlot = '';

% Neurons from which hemisphere to plot. Defined by the somaSide field in
% OF struct array. Type 'Both' to plot all neurons. In that case, responses
% from Right hemisphere will be flipped to match the Left hemisphere neural
% responses and plotted together.
% Accepted arguments: ['L','R','Both'] Default: sideToPlot = 'Both';
sideToPlot = 'Both';

% specify normalization type. Type 'Raw' for no normalization.
% 'Max' will normalize the data per ROI to the maximum response recorded in
% that ROI (see variable maxes).
% 'Zscore' will zscore the data per ROI based on the mean and standard
% deviation of the entire recording session done on that ROI (see variables
% mu and sigma). 'Znorm' normalizes Z-score data by the max Z-score
% observed per fly per neuron
% Accepted arguments: ['Max','Zscore','Raw','Znorm'] Default: 'Raw';
normType = 'Zscore';
%----------------------------end of user input-----------------------------

clearvars l lh h stats

if isempty(stimToPlot)
    stimToPlot = 1:length(Stimuli);
end
if isempty(neuronsToPlot)
    neuronsToPlot = 1:length(regions);
end

prepIndices = cell(length(regions),1);
gcIndices = cell(length(regions),1);
sideIndices = cell(length(regions),1);
prepIndicesFly = cell(length(regions),1);
gcIndicesFly = cell(length(regions),1);
sideIndicesFly = cell(length(regions),1);

if ~isempty(prepToPlot)
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        if isempty(prepIndices{r})
            prepIndices{r} = ...
                contains(PrepTypePerROI{r},prepToPlot);
            prepIndicesFly{r} = ...
                contains(PrepType{r},prepToPlot);
        else
            prepIndices{r} = (prepIndices{r} | ...
                contains(PrepTypePerROI{r},prepToPlot));
            prepIndicesFly{r} = (prepIndicesFly{r} | ...
                contains(PrepType{r},prepToPlot));
        end
    end
else
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        prepIndices{r} = true(size(PrepTypePerROI{r}));
        prepIndicesFly{r} = true(size(PrepType{r}));
    end
end
if ~isempty(gcampsToPlot)
    for p = 1:length(gcampsToPlot) %for each selected GCaMP type
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            if isempty(gcIndices{r})
                gcIndices{r} = ...
                    ismember(GCaMPTypePerROI{r},gcampsToPlot{p});
                gcIndicesFly{r} = ...
                    ismember(GCaMPType{r},gcampsToPlot{p});
            else
                gcIndices{r} = (gcIndices{r} | ...
                    ismember(GCaMPTypePerROI{r},gcampsToPlot{p}));
                gcIndicesFly{r} = (gcIndicesFly{r} | ...
                    ismember(GCaMPType{r},gcampsToPlot{p}));
            end
        end
    end
else
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        gcIndices{r} = true(size(GCaMPTypePerROI{r}));
        gcIndicesFly{r} = true(size(GCaMPType{r}));
    end
end
switch sideToPlot
    case 'Both'
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            sideIndices{r} = true(size(SomaSidePerROI{r}));
            sideIndicesFly{r} = true(size(SomaSide{r}));
        end
        switch normType
            case 'Zscore'
                aucToPlot = Auc_1st_bs_Z_R;
                dataToPlot = OFRes_1st_Z_R;
            case 'Max'
                dataToPlot = Auc_1st_normM;
            case 'Raw'
                dataToPlot1st = Auc_1st_bs;
                dataToPlot2nd = Auc_2nd_bs;
            otherwise
                error(['Normalization type not specified correctly. ',...
                    'Use normType=[''Max'',''Zscore'', ',...
                    'or ''Raw'']'])
        end
    case {'L','R'}
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            if isempty(sideIndices{r})
                sideIndices{r} = ...
                    contains(SomaSidePerROI{r},sideToPlot);
                sideIndicesFly{r} = ...
                    contains(SomaSide{r},sideToPlot);
            else
                sideIndices{r} = (sideIndices{r} | ...
                    contains(SomaSidePerROI{r},sideToPlot));
                sideIndicesFly{r} = (sideIndicesFly{r} | ...
                    contains(SomaSide{r},sideToPlot));
            end
        end
        switch normType
            case 'Zscore'
                dataToPlot1st = Auc_1st_bs_Z;
                dataToPlot2nd = Auc_2nd_bs_Z;
            case 'Max'
                dataToPlot1st = Auc_1st_normM;
                dataToPlot2nd = Auc_2nd_normM;
            case 'Raw'
                dataToPlot1st = Auc_1st_bs;
                dataToPlot2nd = Auc_2nd_bs;
            otherwise
                error(['Normalization type not specified correctly. ',...
                    'Use normType=[''Max'',''Zscore'', ',...
                    'or ''Raw'']'])
        end
    otherwise
        error(['Plotting side not specified correctly. ',...
            'Use sideToPlot = [''L'',''R'' or ''Both'']'])
end


StimXLabel = {'zero','low','int.','fast'};
limsY = [-0.2 1.2];
sp=0.5;
f=figure('WindowState','maximized');
x = zeros(size(neuronsToPlot,2),size(stimToPlotCW,2));
DSI = cell(size(neuronsToPlot,2),size(stimToPlotCW,2));
DSR = cell(size(neuronsToPlot,2),size(stimToPlotCW,2));
DSRnormM = cell(size(neuronsToPlot,2),size(stimToPlotCW,1));
DSRnormMA = cell(size(neuronsToPlot,2),size(stimToPlotCW,1));
my = zeros(1,size(neuronsToPlot,2));

for n = 1:length(neuronsToPlot)
    s = neuronsToPlot(n);
    iPlot = prepIndices{s} & gcIndices{s} & sideIndices{s};
    for k = 1:length(stimToPlotCW)
        mCW = stimToPlotCW(k);
        mCCW = stimToPlotCCW(k);       
        switch regions{s}
            case {'H2axon'}
                rNPD = aucToPlot{s,mCCW}(:,iPlot);
                rPD = aucToPlot{s,mCW}(:,iPlot);
            otherwise
                rPD = aucToPlot{s,mCW}(:,iPlot);
                rNPD = aucToPlot{s,mCCW}(:,iPlot);
        end
        DSI{n,k}=(rPD-rNPD)./(abs(rPD)+abs(rNPD));
        DSR{n,k}=(rPD-rNPD);
        % save the maximum y values for significance plot alignment
        maxbl = nanmax(DSR{n,k});
        if isempty(maxbl)
            maxbl = 0;
        end
        my(n) = max(my(n),maxbl);
    end
    maxnorm = zeros(1,size(DSR{n,k},2));
    for i = 1:size(DSR{n,k},2)
        for k = 1:length(stimToPlotCW)
            if i <= size(DSR{n,k},2)
                maxnorm(1,i) = max([maxnorm(1,i), abs(DSR{n,k}(i))]);
            end
        end
    end
    for k = 1:length(stimToPlotCW)
        DSRnormM{n,k}=DSR{n,k}./maxnorm(1,1:size(DSR{n,k},2));
        DSRnormMA{n,k}=abs(DSR{n,k})./maxnorm(1,1:size(DSR{n,k},2));
        DSI{n,k}(isnan(DSI{n,k})) = [];
        DSR{n,k}(isnan(DSR{n,k})) = [];
        DSRnormM{n,k}(isnan(DSRnormM{n,k})) = [];
        DSRnormMA{n,k}(isnan(DSRnormMA{n,k})) = [];
    end
    for k = 1:4
        sti = k+5;
        x(n,k) = (n-1)*length(stimToPlot) + 1 + (k-1)*sp;
        h(n,k) = notBoxPlot(DSRnormMA{n,sti},x(n,k));hold on
        set(h(n,k).data,'MarkerFaceColor',cReg(s,:))
    end
    % Perform Mann-Whitney-Wilcoxon non parametric test and plot
    % significance indicators if p<0.05
    nck = nchoosek(1:k,2);
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
        %         stats{ii} = mwwtest(DSI{n,nck(ii,1)},DSI{n,nck(ii,2)},0);
        stats{ii} = mwwtest(DSRnormMA{n,nck(ii,1)},DSRnormMA{n,nck(ii,2)},0);
        % uncomment below if you want to display the p-values per pair
        %         disp(['n=',num2str(n),'-',regions{neuronsToPlot(nck(ii,1))},'-',...
        %             regions{neuronsToPlot(nck(ii,2))},'-',num2str(stats{ii}.p(2))])
        if stats{ii}.p < 0.05
            %                     sigline([x(n,nck(ii,1)),x(n,nck(ii,2))],[],my(n)+my(n)*0.03*ii);
        end
    end
end
limsX = xlim;
line(limsX,[0,0],'LineStyle','--','Color','k','LineWidth',1)
set(gca,'XTick',unique(x),'XTickLabel',StimXLabel,...
    'FontSize',14,'TickLabelInterpreter','none',...
    'XTickLabelRotation',45)
set(gca,'YLim',limsY)
ylabel('Asymmetry Index')

for m=1:length(neuronsToPlot)
    hold on
    l(m)=plot(nan,nan,'Color',cReg(neuronsToPlot(m),:),'LineWidth',4);
end
set(gca,'XLim',limsX)
lh=legend(l,regions(neuronsToPlot),'Box','off','FontSize',20,...
    'Location','Best','Interpreter','none');


%% Extended Data Figure 1a
% --------------------------- User Input ----------------------------------
% which OF stimuli to plot. Defined by their index in the Stimuli cell
% array. Plotted in the same order. Leave blank to plot all stimuli.
stimToPlot = find(contains(StimuliSTPooled,'DotYaw_4pxsq_CCW'));

% which neurons to plot. Defined by their index in the regions cell array.
% Plotted in the same order. Leave blank to plot all neurons.
neuronsToPlot = [8,3]; % HS DNp15

% Which calcium indicator(s) to plot. Defined by the genotype field in OF
% struct array. Leave blank to plot everything. Note that it is possible to
% pool data from multiple indicators (e.g. gcampsToPlot = {GC6f','GC7f'})
% Accepted arguments: {'GC6f','GC7f','sytGC7f'} Default: gcampsToPlot = {};
gcampsToPlot = {'GC7f','sytGC7f'};

% Which prep types to plot. Defined by the prepType field in OF struct
% array. Leave blank to plot all experiments.
% Accepted arguments: ['walking','immobilized'] Default: prepToPlot = '';
prepToPlot = 'immobilized';

% Neurons from which hemisphere to plot. Defined by the somaSide field in
% OF struct array. Type 'Both' to plot all neurons. In that case, responses
% from Right hemisphere will be flipped to match the Left hemisphere neural
% responses and plotted together.
% Accepted arguments: ['L','R','Both'] Default: sideToPlot = 'Both';
sideToPlot = 'Both';

% specify normalization type. Type 'Raw' for no normalization.
% 'Max' will normalize the data per ROI to the maximum response recorded in
% that ROI (see variable maxes).
% 'Zscore' will zscore the data per ROI based on the mean and standard
% deviation of the entire recording session done on that ROI (see variables
% mu and sigma). 'Znorm' normalizes Z-score data by the max Z-score
% observed per fly per neuron
% Accepted arguments: ['Max','Zscore','Raw','Znorm'] Default: 'Raw';
normType = 'Zscore';

% Specify plot type. Type 'Line' for plotting each fly with a different
% colored line together with the mean trace. 'Shaded' will plot the std as
% a shaded area around the mean trace instead.
% Accepted arguments: ['Line','Shaded']
plotType = 'Shaded';

switch normType
    case 'Max'
        limsY = [-0.5 1.1]; % Y-axis limits
    case 'Raw'
        limsY = [-3 10]; % Y-axis limits
    case 'Zscore'
        switch plotType
            case 'Line'
                limsY = [-5 30]; % Y-axis limits
            case 'Shaded'
                %                 limsY = [-3.5 30]; % Y-axis limits
                %                 limsY = [-5 15]; % Y-axis limits
                limsY = [-5 35]; % Y-axis limits
        end
    case 'Znorm'
        limsY = [-0.25 1]; % Y-axis limits
end
%----------------------------end of user input-----------------------------

% Define the subset of data that will be plotted. This will be used in the
% plotSelectDFF function to plot only the specified subset of data as well
% as to calculate and display number of ROIs being plotted
if isempty(stimToPlot)
    stimToPlot = 1:length(StimuliSTPooled);
end
if isempty(neuronsToPlot)
    neuronsToPlot = 1:length(regions);
end
prepIndices = cell(length(regions),1);
gcIndices = cell(length(regions),1);
sideIndices = cell(length(regions),1);
prepIndicesFly = cell(length(regions),1);
gcIndicesFly = cell(length(regions),1);
sideIndicesFly = cell(length(regions),1);
if ~isempty(prepToPlot)
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        if isempty(prepIndices{r})
            prepIndices{r} = ...
                contains(PrepTypePerROI{r},prepToPlot);
            prepIndicesFly{r} = ...
                contains(PrepType{r},prepToPlot);
        else
            prepIndices{r} = (prepIndices{r} | ...
                contains(PrepTypePerROI{r},prepToPlot));
            prepIndicesFly{r} = (prepIndicesFly{r} | ...
                contains(PrepType{r},prepToPlot));
        end
    end
else
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        prepIndices{r} = true(size(PrepTypePerROI{r}));
        prepIndicesFly{r} = true(size(PrepType{r}));
    end
end
if ~isempty(gcampsToPlot)
    for p = 1:length(gcampsToPlot) %for each selected GCaMP type
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            if isempty(gcIndices{r})
                gcIndices{r} = ...
                    ismember(GCaMPTypePerROI{r},gcampsToPlot{p});
                gcIndicesFly{r} = ...
                    ismember(GCaMPType{r},gcampsToPlot{p});
            else
                gcIndices{r} = (gcIndices{r} | ...
                    ismember(GCaMPTypePerROI{r},gcampsToPlot{p}));
                gcIndicesFly{r} = (gcIndicesFly{r} | ...
                    ismember(GCaMPType{r},gcampsToPlot{p}));
            end
        end
    end
else
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        gcIndices{r} = true(size(GCaMPTypePerROI{r}));
        gcIndicesFly{r} = true(size(GCaMPType{r}));
    end
end
switch sideToPlot
    case 'Both'
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            sideIndices{r} = true(size(SomaSidePerROI{r}));
            sideIndicesFly{r} = true(size(SomaSide{r}));
        end
        switch normType
            case 'Zscore'
                dataToPlot = STResP_Z;
            case 'Max'
                dataToPlot = STResP_normM;
            case 'Raw'
                dataToPlot = STResP;
            case 'Znorm'
                %                 dataToPlot1st = OFRes_1st_Znorm;
                %                 dataToPlot2nd = OFRes_2nd_Znorm;
            otherwise
                error(['Normalization type not specified correctly. ',...
                    'Use normType=[''Max'',''Zscore'',''Znorm'' ',...
                    'or ''Raw'']'])
        end
    case {'L','R'}
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            if isempty(sideIndices{r})
                sideIndices{r} = ...
                    contains(SomaSidePerROI{r},sideToPlot);
                sideIndicesFly{r} = ...
                    contains(SomaSide{r},sideToPlot);
            else
                sideIndices{r} = (sideIndices{r} | ...
                    contains(SomaSidePerROI{r},sideToPlot));
                sideIndicesFly{r} = (sideIndicesFly{r} | ...
                    contains(SomaSide{r},sideToPlot));
            end
        end
        switch normType
            case 'Zscore'
                dataToPlot = STResP_Raw_Z;
            case 'Max'
                dataToPlot = STResP_Raw_normM;
            case 'Raw'
                dataToPlot1st = STResP_Raw;
            case 'Znorm'
                %                 dataToPlot1st = OFRes_1st_Znorm;
                %                 dataToPlot2nd = OFRes_2nd_Znorm;
            otherwise
                error(['Normalization type not specified correctly. ',...
                    'Use normType=[''Max'',''Zscore'',''Znorm'' ',...
                    'or ''Raw'']'])
        end
    otherwise
        error(['Plotting side not specified correctly. ',...
            'Use sideToPlot = [''L'',''R'' or ''Both'']'])
end
% cReg = [0 0 1;1 0 1;0.4 0 0.4;1 0 0;0.4 0 0;0 1 0;0 1 1;1 0.5 0];
f=plotSelectDFFST(dataToPlot,neuronsToPlot,...
    stimToPlot,minFrameRate,cReg,prepIndices,prepIndicesFly,...
    gcIndices,gcIndicesFly,sideIndices,sideIndicesFly,StimuliSTPooled,...
    regions,limsY,plotType);
%% Extended Data Figure 1b
neuronsToPlot = [3,8];
PlotType = 'bar'; % accepted arguments  = {'bar','shaded'}

% initialize speed tuning curve related variables
STcurves = cell(length(regions),1);
STcurves_Z = cell(length(regions),1);
STslowfast = cell(length(regions),2);

% create a new figure
f=figure('WindowState','maximized');
hold on
for i =  neuronsToPlot % for each unique neuron-region combination
    switch regions{i} % define which OF ST curve will be plotted per neuron
        case {'H2axon'}
            stimIndex = find(contains(StimuliSTPooled,'DotYaw_4pxsq_CW'));
        otherwise
            stimIndex = find(contains(StimuliSTPooled,'DotYaw_4pxsq_CCW'));
    end
    stimNames = StimuliSTPooled(stimIndex);

    % initialize ST curve related variables
    speeds = cell(1,length(stimNames));
    STcurves{i} = nan(size(STResP_Z{i,stimIndex(1)},2),length(stimNames));

    % Define the indices for each trial. Used to define the range for auc
    % calculations (seconds 1-2 for baseline, 2-6 for stimulus-evoked auc)
    frameNum = size(STResP_Z{i,1},1);
    TimeIndices = linspace(0,frameNum/(minFrameRate/1000),frameNum);

    for spd=1:length(stimNames) % for each stimulus speed
        % get the value of speeds used for the x-axis label
        dashes = strfind(stimNames{spd},'_');
        dashes2 = strfind(stimNames{spd},'px/s');
        tempspeed = stimNames{spd}(dashes(end)+1:dashes2(end)-1);
        speeds{spd} = str2double(tempspeed)*2.25;
        if ~isempty(STResP_Z{i,stimIndex(spd)}) % if there is imaging data
            if ~isempty(STResP_Z{i,stimIndex(spd)}...
                    (~isnan(STResP_Z{i,stimIndex(spd)})))%if data has value
                % auc of the stimulus evoked calcium response (seconds 2-6)
                STcurves{i}(:,spd) = trapz(STResP_Z{i,stimIndex(spd)}...
                    ((find(TimeIndices<2,1,'last')+1):end,:));
                STcurves{i}(:,spd) = trapz(STResP_Z{i,stimIndex(spd)}...
                    ((find(TimeIndices<2,1,'last')+1):...
                    (find(TimeIndices<4,1,'last')),:));
            end
        end
    end
    STcurves{i}(isnan(STcurves{i}(:,spd)),:) = []; % remove NaN data points
    % normalize each ROI based on its absolute maximum auc response
    STcurves_Z{i} = STcurves{i}./max(abs(STcurves{i}),[],2);
    switch PlotType
        case 'shaded'
            % plot the curves with shaded error bars
            shadedErrorBar(1:length(speeds),mean(STcurves_Z{i}),...
                std(STcurves_Z{i})/sqrt(size(STcurves_Z{i},1)),'lineprops',...
                {'Color',cReg(i,:),'LineWidth',3},'transparent',1)
        case 'bar'
            % plot the curves with line error bars
            eh = errorbar(1:length(speeds),mean(STcurves_Z{i}),...
                std(STcurves_Z{i})/sqrt(size(STcurves_Z{i},1)));
            set(eh, 'Color', cReg(i,:),'LineWidth',3, 'CapSize', 15)
        otherwise
            error('unrecognized PlotType variable')
    end
end

limsX = [0 length(speeds)+1];
% modify X-axis labels to match the speeds presented, add y-label
set(gca,'XTick',1:length(speeds),'XTickLabel',speeds,'XLim',...
    limsX,'FontSize',18)
xlabel('Stimulus speed (deg/s)');
ylabel('Normalized response');
% add a dashed line at zero for reference
line(limsX,[0,0],'LineStyle','--','Color','k','LineWidth',1)
% add imaginary lines for each neuron (used to make the legend)
l = gobjects(1,length(neuronsToPlot));
for k=1:length(neuronsToPlot)
    hold on
    l(k)=plot(nan,nan,'Color',cReg(neuronsToPlot(k),:),'LineWidth',4);
end
lh=legend(l,regions(neuronsToPlot),'Box','off','FontSize',24); % add legend
clearvars l lh eh

%% Extended Data Figure 1c
% --------------------------- User Input ----------------------------------
% which OF stimuli to plot. Defined by their index in the Stimuli cell
% array. Plotted in the same order. Leave blank to plot all stimuli.
stimToPlot = [2 10 15 17 16 14];

% which neurons to plot. Defined by their index in the regions cell array.
% Plotted in the same order. Leave blank to plot all neurons.
neuronsToPlot = [8,3];

% Which calcium indicator(s) to plot. Defined by the genotype field in OF
% struct array. Leave blank to plot everything. Note that it is possible to
% pool data from multiple indicators (e.g. gcampsToPlot = {GC6f','GC7f'})
% Accepted arguments: {'GC6f','GC7f','sytGC7f','GC8f'}
% Default: gcampsToPlot = {};
gcampsToPlot = {'GC7f','sytGC7f'};

% Which prep types to plot. Defined by the prepType field in OF struct
% array. Leave blank to plot all experiments. Keep in mind that all of the
% OF sensitivity experiments in walking prep were done when the fly was
% stationary and not moving/grooming/struggling etc.
% Accepted arguments: ['walking','immobilized'] Default: prepToPlot = '';
prepToPlot = 'immobilized';

% Neurons from which hemisphere to plot. Defined by the somaSide field in
% OF struct array. Type 'Both' to plot all neurons. In that case, responses
% from Right hemisphere will be flipped to match the Left hemisphere neural
% responses and plotted together.
% Accepted arguments: ['L','R','Both'] Default: sideToPlot = 'Both';
sideToPlot = 'Both';

% specify normalization type. Type 'DFF' for (F-Fo)/Fo normalization.
% 'Max' will normalize the data per ROI to the maximum DFF response
% recorded in that ROI (see variable maxes).
% 'Zscore' will z-score the data per ROI based on the baseline mean and
% standard deviation for each trial done on that ROI.
% 'Znorm' normalizes Z-score data by the max Z-score recorded from each ROI
% Accepted arguments: ['DFF','Max','Zscore','Znorm']
normType = 'Zscore';

% Specify plot type. Type 'Line' for plotting each fly with a different
% colored line together with the mean trace. 'Shaded' will plot the std as
% a shaded area around the mean trace instead.
% Accepted arguments: ['Line','Shaded']
plotType = 'Shaded';

% Define Y-axis limits for each normalization type.
switch normType
    case 'Max'
        limsY = [-0.5 1.1]; % Y-axis limits
    case 'DFF'
        limsY = [-3 10]; % Y-axis limits
    case 'Zscore'
        switch plotType
            case 'Line'
                limsY = [-6 50]; % Y-axis limits
            case 'Shaded'
                limsY = [-3.5 28]; % Y-axis limits
        end
    case 'Znorm'
        limsY = [-0.4 1.05]; % Y-axis limits
end
%----------------------------end of user input-----------------------------
H2flip = 1;
H2i = find(contains(regions,'H2axon')); % index of H2
% Define the subset of data that will be plotted. This will be used in the
% plotSelectDFF function to plot only the specified subset of data as well
% as to calculate and display number of flies & ROIs being plotted
if isempty(stimToPlot)
    stimToPlot = 1:length(Stimuli);
end
if isempty(neuronsToPlot)
    neuronsToPlot = 1:length(regions);
end
prepIndices = cell(length(regions),1);
gcIndices = cell(length(regions),1);
sideIndices = cell(length(regions),1);
prepIndicesFly = cell(length(regions),1);
gcIndicesFly = cell(length(regions),1);
sideIndicesFly = cell(length(regions),1);
if ~isempty(prepToPlot)
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        if isempty(prepIndices{r})
            prepIndices{r} = ...
                contains(PrepTypePerROI{r},prepToPlot);
            prepIndicesFly{r} = ...
                contains(PrepType{r},prepToPlot);
        else
            prepIndices{r} = (prepIndices{r} | ...
                contains(PrepTypePerROI{r},prepToPlot));
            prepIndicesFly{r} = (prepIndicesFly{r} | ...
                contains(PrepType{r},prepToPlot));
        end
    end
else
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        prepIndices{r} = true(size(PrepTypePerROI{r}));
        prepIndicesFly{r} = true(size(PrepType{r}));
    end
end
if ~isempty(gcampsToPlot)
    for p = 1:length(gcampsToPlot) %for each selected GCaMP type
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            if isempty(gcIndices{r})
                gcIndices{r} = ...
                    ismember(GCaMPTypePerROI{r},gcampsToPlot{p});
                gcIndicesFly{r} = ...
                    ismember(GCaMPType{r},gcampsToPlot{p});
            else
                gcIndices{r} = (gcIndices{r} | ...
                    ismember(GCaMPTypePerROI{r},gcampsToPlot{p}));
                gcIndicesFly{r} = (gcIndicesFly{r} | ...
                    ismember(GCaMPType{r},gcampsToPlot{p}));
            end
        end
    end
else
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        gcIndices{r} = true(size(GCaMPTypePerROI{r}));
        gcIndicesFly{r} = true(size(GCaMPType{r}));
    end
end
switch sideToPlot
    case 'Both'
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            sideIndices{r} = true(size(SomaSidePerROI{r}));
            sideIndicesFly{r} = true(size(SomaSide{r}));
        end
        switch normType
            case 'Zscore'
                dataToPlot1st = OFRes_1st_Z_R;
                dataToPlot2nd = OFRes_2nd_Z_R;
                if H2flip
                    dataToPlot1st(H2i,:) = OFRes_1st_Z(H2i,:);
                    dataToPlot2nd(H2i,:) = OFRes_2nd_Z(H2i,:);
                end
            case 'Max'
                dataToPlot1st = OFRes_1st_normM;
                dataToPlot2nd = OFRes_2nd_normM;
            case 'DFF'
                dataToPlot1st = OFRes_1st;
                dataToPlot2nd = OFRes_2nd;
            case 'Znorm'
                dataToPlot1st = OFRes_1st_Znorm;
                dataToPlot2nd = OFRes_2nd_Znorm;
            otherwise
                error(['Normalization type not specified correctly. ',...
                    'Use normType=[''Max'',''Zscore'',''Znorm'' ',...
                    'or ''Raw'']'])
        end
    case {'L','R'}
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            if isempty(sideIndices{r})
                sideIndices{r} = ...
                    contains(SomaSidePerROI{r},sideToPlot);
                sideIndicesFly{r} = ...
                    contains(SomaSide{r},sideToPlot);
            else
                sideIndices{r} = (sideIndices{r} | ...
                    contains(SomaSidePerROI{r},sideToPlot));
                sideIndicesFly{r} = (sideIndicesFly{r} | ...
                    contains(SomaSide{r},sideToPlot));
            end
        end
        switch normType
            case 'Zscore'
                dataToPlot1st = OFRes_Raw_1st_Z;
                dataToPlot2nd = OFRes_Raw_2nd_Z;
            case 'Max'
                dataToPlot1st = OFRes_Raw_1st_normM;
                dataToPlot2nd = OFRes_Raw_2nd_normM;
            case 'Raw'
                dataToPlot1st = OFRes_Raw_1st;
                dataToPlot2nd = OFRes_Raw_2nd;
            case 'Znorm'
                dataToPlot1st = OFRes_1st_Znorm;
                dataToPlot2nd = OFRes_2nd_Znorm;
            otherwise
                error(['Normalization type not specified correctly. ',...
                    'Use normType=[''Max'',''Zscore'',''Znorm'' ',...
                    'or ''Raw'']'])
        end
    otherwise
        error(['Plotting side not specified correctly. ',...
            'Use sideToPlot = [''L'',''R'' or ''Both'']'])
end

f=plotSelectDFF(dataToPlot1st,dataToPlot2nd,neuronsToPlot,...
    stimToPlot,minFrameRate,cReg,prepIndices,prepIndicesFly,...
    gcIndices,gcIndicesFly,sideIndices,sideIndicesFly,Stimuli,...
    regions,limsY,plotType,normType);

%% Extended Data Figure 1d

% --------------------------- User Input ----------------------------------
neuronsToPlot = [8,3];
% Specify plot type. Type 'Line' for plotting each fly with a different
% colored line together with the mean trace. 'Shaded' will plot the std as
% a shaded area around the mean trace instead.
% Accepted arguments: ['Line','Shaded']
plotType = 'Shaded';
%----------------------------end of user input-----------------------------

% Translational stimuli used (as of 210421)
StimTR = {'Sideslip','ProgSide315','Progressive','ProgSide45'};

% Find the corresponding indices for the stimuli defined in StimTR in the
% Stimuli cell array
iTR = nan(1,length(StimTR));
for i = 1:length(iTR)
    iTR(i)=find(contains(Stimuli,StimTR(i))&~contains(Stimuli,'Left') & ...
        ~contains(Stimuli,'Right'));
end

lmY = [-3.5 30]; % Y-axis limits for HS
cRect = [0 0.5 0 0.1]; % color of the rectangle depicting stimulus motion
f=figure('WindowState','maximized'); hold on
for i = 1:length(neuronsToPlot) % for each unique neuron-region combination
    n = neuronsToPlot(i);
    frameNum = size(OFRes_1st_Z_R{n,1},1);
    TimeAxisPD = linspace(0,frameNum/(minFrameRate/1000),frameNum);
    %     limsX = [min(TimeAxisPD) ceil(max(TimeAxisPD))];
    limsX = [min(TimeAxisPD)+1 ceil(max(TimeAxisPD))-0.5];

    for j = 1:length(StimTR) % for each OF stimulus
        if ~isempty(OFRes_1st_Z_R{n,iTR(j)}) % if there is imaging data
            CounterTRn = sum(~isnan(OFRes_1st_Z_R{n,iTR(j)}(1,:)));
            switch StimTR{j}
                case 'Sideslip'
                    subplot(3,3,6)
                case 'ProgSide315'
                    subplot(3,3,1)
                case 'Progressive'
                    subplot(3,3,2)
                case 'ProgSide45'
                    subplot(3,3,3)
            end
            % plot individual neurites, and mean across neurites, 1st epoch
            switch plotType
                case 'Line'
                    plot(TimeAxisPD,OFRes_1st_Z_R{n,iTR(j)}); hold on;
                    plot(TimeAxisPD,nanmean(OFRes_1st_Z_R{n,iTR(j)},2),...
                        'Color',cReg(n,:),'LineWidth',4);
                case 'Shaded'
                    h=shadedErrorBar(TimeAxisPD,...
                        nanmean(OFRes_1st_Z_R{n,iTR(j)},2),...
                        nanstd(OFRes_1st_Z_R{n,iTR(j)},0,2)/...
                        sqrt(CounterTRn));hold on
                    set(h.mainLine,'Color',cReg(n,:),...
                        'LineWidth',4)
                    set(h.patch,'FaceColor',cReg(n,:),...
                        'FaceAlpha',0.2)
                    set(h.edge(1),'Color',cReg(n,:))
                    set(h.edge(2),'Color',cReg(n,:))
                    clearvars h
            end

            rectangle('Position',[2 lmY(1) 2 lmY(2)-lmY(1)],...
                'FaceColor', cRect);
            set(gca,'YLim',lmY,'XLim',limsX)
            switch StimTR{j}
                case 'ProgSide315'
                    ylabel('Response (z-score)')
                    xlabel('Time (sec)')
                    set(gca,'FontSize',18)
                otherwise
                    axis off
            end
            line(limsX,[0 0],'Color','k','LineWidth',0.5,'LineStyle','--')
        end

        if  ~isempty(OFRes_2nd_Z_R{n,iTR(j)}) % if there is imaging data
            switch StimTR{j}
                case 'Sideslip'
                    subplot(3,3,4)
                case 'ProgSide315'
                    subplot(3,3,9)
                case 'Progressive'
                    subplot(3,3,8)
                case 'ProgSide45'
                    subplot(3,3,7)
            end
            % plot individual neurites, and mean across neurites, 2nd epoch
            switch plotType
                case 'Line'
                    plot(TimeAxisPD,OFRes_2nd_Z_R{n,iTR(j)}); hold on;
                    plot(TimeAxisPD,nanmean(OFRes_2nd_Z_R{n,iTR(j)},2),...
                        'color',cReg(n,:),'LineWidth',4);
                case 'Shaded'
                    h=shadedErrorBar(TimeAxisPD,...
                        nanmean(OFRes_2nd_Z_R{n,iTR(j)},2),...
                        nanstd(OFRes_2nd_Z_R{n,iTR(j)},0,2)/...
                        sqrt(CounterTRn));hold on
                    set(h.mainLine,'Color',cReg(n,:),...
                        'LineWidth',4)
                    set(h.patch,'FaceColor',cReg(n,:),...
                        'FaceAlpha',0.2)
                    set(h.edge(1),'Color',cReg(n,:))
                    set(h.edge(2),'Color',cReg(n,:))
                    clearvars h
            end
            rectangle('Position',[2 lmY(1) 2 lmY(2)-lmY(1)],...
                'FaceColor', cRect);
            set(gca,'YLim',lmY,'XLim',limsX)
            axis off
            line(limsX,[0 0],'Color','k','LineWidth',0.5,'LineStyle','--')
        end
    end
end
%% Extended Data Figure 1e

% --------------------------- User Input ----------------------------------
% which neurons to plot. Defined by their index in the regions cell array.
% Plotted in the same order. Leave blank to plot all neurons.

neuronsToPlot = [8,3];

% Which calcium indicator(s) to plot. Defined by the genotype field in OF
% struct array. Leave blank to plot everything. Note that it is possible to
% pool data from multiple indicators (e.g. gcampsToPlot = {GC6f','GC7f'})
% Accepted arguments: {'GC6f','GC7f','sytGC7f'} Default: gcampsToPlot = {};
gcampsToPlot = {'sytGC7f','GC7f'};

% Which prep types to plot. Defined by the prepType field in OF struct
% array. Leave blank to plot all experiments.
% Accepted arguments: ['walking','immobilized'] Default: prepToPlot = '';
prepToPlot = 'immobilized';

% Neurons from which hemisphere to plot. Defined by the somaSide field in
% OF struct array. Type 'Both' to plot all neurons. In that case, responses
% from Right hemisphere will be flipped to match the Left hemisphere neural
% responses and plotted together.
% Accepted arguments: ['L','R','Both'] Default: sideToPlot = 'Both';
sideToPlot = 'Both';

% specify normalization type. Type 'Raw' for no normalization.
% 'Max' will normalize the data per ROI to the maximum response recorded in
% that ROI (see variable maxes).
% 'Zscore' will zscore the data per ROI based on the mean and standard
% deviation of the entire recording session done on that ROI (see variables
% mu and sigma). 'Znorm' normalizes Z-score data by the max Z-score
% observed per fly per neuron
% Accepted arguments: ['Zscore','Raw'] Default: 'Raw';
normType = 'Zscore';

% Specify plot type. Type 'Box' for plotting each fly with a different
% colored dot together with the population mean and sd as boxes.
% Currently only 'Box' is supported.
% Accepted arguments: ['Box']
plotType = 'Box';

% Specify grouping. Type 'neuron' for plotting all stimuli for each neuron
% grouped together. 'stimulus' will group all neurons for each stimulus
% Accepted arguments: ['neuron','stimulus']
groupType = 'neuron';

%----------------------------end of user input-----------------------------

peakbin = 500; % time bin for the averaged peak response (ms)
H2Flip = 0;
iH2 = find(contains(regions,'H2axon')); % index of H2
% index for recordings that will come from the right hemisphere
[~,loc] = ismember(StimuliFlipped,Stimuli);

clearvars l lh h stats
stm = {'Yaw_','Progressive','YawRightFlipped','RollInv','Pitch'};
auci = find(contains(Stimuli,stm) & BinocularIndices);
yi = auci(contains(Stimuli(auci),'Yaw_'));
yfi = auci(contains(Stimuli(auci),'YawRightFlipped'));
pri = auci(contains(Stimuli(auci),'Progressive'));
ssi = find(contains(Stimuli,'Sideslip'));
ps45i = find(contains(Stimuli,'ProgSide45'));

% Define the subset of data that will be plotted. This will be used in the
% plotSelectAUC function to plot only the specified subset of data as well
% as to calculate and display number of flies & ROIs being plotted

if isempty(neuronsToPlot)
    neuronsToPlot = 1:length(regions);
end

prepIndices = cell(length(regions),1);
gcIndices = cell(length(regions),1);
sideIndices = cell(length(regions),1);
prepIndicesFly = cell(length(regions),1);
gcIndicesFly = cell(length(regions),1);
sideIndicesFly = cell(length(regions),1);

if ~isempty(prepToPlot)
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        if isempty(prepIndices{r})
            prepIndices{r} = ...
                contains(PrepTypePerROI{r},prepToPlot);
            prepIndicesFly{r} = ...
                contains(PrepType{r},prepToPlot);
        else
            prepIndices{r} = (prepIndices{r} | ...
                contains(PrepTypePerROI{r},prepToPlot));
            prepIndicesFly{r} = (prepIndicesFly{r} | ...
                contains(PrepType{r},prepToPlot));
        end
    end
else
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        prepIndices{r} = true(size(PrepTypePerROI{r}));
        prepIndicesFly{r} = true(size(PrepType{r}));
    end
end
if ~isempty(gcampsToPlot)
    for p = 1:length(gcampsToPlot) %for each selected GCaMP type
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            if isempty(gcIndices{r})
                gcIndices{r} = ...
                    ismember(GCaMPTypePerROI{r},gcampsToPlot{p});
                gcIndicesFly{r} = ...
                    ismember(GCaMPType{r},gcampsToPlot{p});
            else
                gcIndices{r} = (gcIndices{r} | ...
                    ismember(GCaMPTypePerROI{r},gcampsToPlot{p}));
                gcIndicesFly{r} = (gcIndicesFly{r} | ...
                    ismember(GCaMPType{r},gcampsToPlot{p}));
            end
        end
    end
else
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        gcIndices{r} = true(size(GCaMPTypePerROI{r}));
        gcIndicesFly{r} = true(size(GCaMPType{r}));
    end
end
switch sideToPlot
    case 'Both'
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            sideIndices{r} = true(size(SomaSidePerROI{r}));
            sideIndicesFly{r} = true(size(SomaSide{r}));
        end
        switch normType
            case 'Zscore'
                dataToPlot1st = OFRes_1st_Z;
                dataToPlot2nd = OFRes_2nd_Z;
            case 'Max'
                dataToPlot1st = Auc_1st_normM;
                dataToPlot2nd = Auc_2nd_normM;
            case 'Raw'
                dataToPlot1st = Auc_1st_bs;
                dataToPlot2nd = Auc_2nd_bs;
            otherwise
                error(['Normalization type not specified correctly. ',...
                    'Use normType=[''Max'',''Zscore'', ',...
                    'or ''Raw'']'])
        end
    case {'L','R'}
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            if isempty(sideIndices{r})
                sideIndices{r} = ...
                    contains(SomaSidePerROI{r},sideToPlot);
                sideIndicesFly{r} = ...
                    contains(SomaSide{r},sideToPlot);
            else
                sideIndices{r} = (sideIndices{r} | ...
                    contains(SomaSidePerROI{r},sideToPlot));
                sideIndicesFly{r} = (sideIndicesFly{r} | ...
                    contains(SomaSide{r},sideToPlot));
            end
        end
        switch normType
            case 'Zscore'
                dataToPlot1st = Auc_1st_bs_Z;
                dataToPlot2nd = Auc_2nd_bs_Z;
            case 'Max'
                dataToPlot1st = Auc_1st_normM;
                dataToPlot2nd = Auc_2nd_normM;
            case 'Raw'
                dataToPlot1st = Auc_1st_bs;
                dataToPlot2nd = Auc_2nd_bs;
            otherwise
                error(['Normalization type not specified correctly. ',...
                    'Use normType=[''Max'',''Zscore'', ',...
                    'or ''Raw'']'])
        end
    otherwise
        error(['Plotting side not specified correctly. ',...
            'Use sideToPlot = [''L'',''R'' or ''Both'']'])
end

frameNum = size(dataToPlot1st{1,1},1);
TimeAxisPD = linspace(0,frameNum/(minFrameRate/1000),frameNum);
iPeak = find(TimeAxisPD>4-peakbin/1000,1,'first'):...
    find(TimeAxisPD<4,1,'last');

compNames = {'Prog-Yaw','Prog-Sideslip','ProgSide45-Sideslip'};
ncomps = length(compNames);
sp=0.5;
x = zeros(ncomps,size(neuronsToPlot,2));
trIndex = cell(ncomps,size(neuronsToPlot,2));
f=figure('WindowState','maximized');
limsY = [-1.2 1.7];
my = zeros(1,ncomps);

for n = 1:ncomps
    for k = neuronsToPlot
        m = find(neuronsToPlot==k);
        x(n,m) = (n-1)*length(neuronsToPlot) + 1 + (m-1)*sp;
        iPlot = prepIndices{k} & gcIndices{k} & sideIndices{k};
        switch compNames{n}
            case 'Prog-Yaw'
                switch regions{k}
                    case {'H2axon'}
                        PKprog = mean(dataToPlot2nd{k,pri}(iPeak,iPlot))-mean(dataToPlot1st{k,pri}(iPeak,iPlot));
                        tPKprog = abs(mean(dataToPlot2nd{k,pri}(iPeak,iPlot))-mean(dataToPlot1st{k,pri}(iPeak,iPlot)));
                        PKyaw = mean(dataToPlot1st{k,yi}(iPeak,iPlot))-mean(dataToPlot2nd{k,yi}(iPeak,iPlot));
                        tPKyaw = abs(mean(dataToPlot1st{k,yi}(iPeak,iPlot))-mean(dataToPlot2nd{k,yi}(iPeak,iPlot)));
                        trIndex{n,m}=(PKprog-PKyaw)./(tPKprog+tPKyaw);
                    otherwise
                        PKprog = mean(dataToPlot1st{k,pri}(iPeak,iPlot))-mean(dataToPlot2nd{k,pri}(iPeak,iPlot));
                        tPKprog = abs(mean(dataToPlot1st{k,pri}(iPeak,iPlot))-mean(dataToPlot2nd{k,pri}(iPeak,iPlot)));
                        PKyaw = mean(dataToPlot2nd{k,yi}(iPeak,iPlot))-mean(dataToPlot1st{k,yi}(iPeak,iPlot));
                        tPKyaw = abs(mean(dataToPlot2nd{k,yi}(iPeak,iPlot))-mean(dataToPlot1st{k,yi}(iPeak,iPlot)));
                        trIndex{n,m}=(PKprog-PKyaw)./(tPKprog+tPKyaw);
                end
            case 'YawFlipped-Yaw'
                switch regions{k}
                    case {'H2axon'}
                        PKyf = mean(dataToPlot1st{k,yfi}(iPeak,iPlot))-mean(dataToPlot2nd{k,yfi}(iPeak,iPlot));
                        tPKyf = abs(mean(dataToPlot2nd{k,yfi}(iPeak,iPlot))-mean(dataToPlot1st{k,yfi}(iPeak,iPlot)));
                        PKyaw = mean(dataToPlot1st{k,yi}(iPeak,iPlot))-mean(dataToPlot2nd{k,yi}(iPeak,iPlot));
                        tPKyaw = abs(mean(dataToPlot1st{k,yi}(iPeak,iPlot))-mean(dataToPlot2nd{k,yi}(iPeak,iPlot)));
                        trIndex{n,m}=(PKyf-PKyaw)./(tPKyf+tPKyaw);
                    case {'H2rnaxon'}
                        PKyf = mean(dataToPlot1st{k,yfi}(iPeak,iPlot))-mean(dataToPlot2nd{k,yfi}(iPeak,iPlot));
                        tPKyf = abs(mean(dataToPlot2nd{k,yfi}(iPeak,iPlot))-mean(dataToPlot1st{k,yfi}(iPeak,iPlot)));
                        PKyaw = mean(dataToPlot2nd{k,yi}(iPeak,iPlot))-mean(dataToPlot1st{k,yi}(iPeak,iPlot));
                        tPKyaw = abs(mean(dataToPlot1st{k,yi}(iPeak,iPlot))-mean(dataToPlot2nd{k,yi}(iPeak,iPlot)));
                        trIndex{n,m}=(PKyf-PKyaw)./(tPKyf+tPKyaw);
                    otherwise
                        PKyf = mean(dataToPlot2nd{k,yfi}(iPeak,iPlot))-mean(dataToPlot1st{k,yfi}(iPeak,iPlot));
                        tPKyf = abs(mean(dataToPlot1st{k,yfi}(iPeak,iPlot))-mean(dataToPlot2nd{k,yfi}(iPeak,iPlot)));
                        PKyaw = mean(dataToPlot2nd{k,yi}(iPeak,iPlot))-mean(dataToPlot1st{k,yi}(iPeak,iPlot));
                        tPKyaw = abs(mean(dataToPlot2nd{k,yi}(iPeak,iPlot))-mean(dataToPlot1st{k,yi}(iPeak,iPlot)));
                        trIndex{n,m}=(PKyf-PKyaw)./(tPKyf+tPKyaw);
                end

            case 'Prog-Sideslip'
                switch regions{k}
                    case {'H2axon'}
                        PKss = mean(dataToPlot2nd{k,ssi}(iPeak,iPlot))-mean(dataToPlot1st{k,ssi}(iPeak,iPlot));
                        tPKss = abs(mean(dataToPlot2nd{k,ssi}(iPeak,iPlot))-mean(dataToPlot1st{k,ssi}(iPeak,iPlot)));
                        PKprog = mean(dataToPlot2nd{k,pri}(iPeak,iPlot))-mean(dataToPlot1st{k,pri}(iPeak,iPlot));
                        tPKprog = abs(mean(dataToPlot2nd{k,pri}(iPeak,iPlot))-mean(dataToPlot1st{k,pri}(iPeak,iPlot)));
                        trIndex{n,m}=(PKprog-PKss)./(tPKprog+tPKss);
                    otherwise
                        PKss = mean(dataToPlot1st{k,ssi}(iPeak,iPlot))-mean(dataToPlot2nd{k,ssi}(iPeak,iPlot));
                        tPKss = abs(mean(dataToPlot1st{k,ssi}(iPeak,iPlot))-mean(dataToPlot2nd{k,ssi}(iPeak,iPlot)));
                        PKprog = mean(dataToPlot1st{k,pri}(iPeak,iPlot))-mean(dataToPlot2nd{k,pri}(iPeak,iPlot));
                        tPKprog = abs(mean(dataToPlot1st{k,pri}(iPeak,iPlot))-mean(dataToPlot2nd{k,pri}(iPeak,iPlot)));
                        trIndex{n,m}=(PKprog-PKss)./(tPKprog+tPKss);
                end

            case 'ProgSide45-Sideslip'
                switch regions{k}
                    case {'H2axon'}
                        PKprog45 = mean(dataToPlot2nd{k,ps45i}(iPeak,iPlot))-mean(dataToPlot1st{k,ps45i}(iPeak,iPlot));
                        tPKprog45 = abs(mean(dataToPlot2nd{k,ps45i}(iPeak,iPlot))-mean(dataToPlot1st{k,ps45i}(iPeak,iPlot)));
                        PKss = mean(dataToPlot2nd{k,ssi}(iPeak,iPlot))-mean(dataToPlot1st{k,ssi}(iPeak,iPlot));
                        tPKss = abs(mean(dataToPlot2nd{k,ssi}(iPeak,iPlot))-mean(dataToPlot1st{k,ssi}(iPeak,iPlot)));
                        trIndex{n,m}=(PKprog45-PKss)./(tPKprog45+tPKss);
                    case {'H2rnaxon'}
                        PKprog45 = mean(dataToPlot1st{k,ps45i}(iPeak,iPlot))-mean(dataToPlot2nd{k,ps45i}(iPeak,iPlot));
                        tPKprog45 = abs(mean(dataToPlot1st{k,ps45i}(iPeak,iPlot))-mean(dataToPlot2nd{k,ps45i}(iPeak,iPlot)));
                        PKss = mean(dataToPlot1st{k,ssi}(iPeak,iPlot))-mean(dataToPlot2nd{k,ssi}(iPeak,iPlot));
                        tPKss = abs(mean(dataToPlot1st{k,ssi}(iPeak,iPlot))-mean(dataToPlot2nd{k,ssi}(iPeak,iPlot)));
                        trIndex{n,m}=(PKprog45-PKss)./(tPKprog45+tPKss);
                    otherwise
                        PKprog45 = mean(dataToPlot1st{k,ps45i}(iPeak,iPlot))-mean(dataToPlot2nd{k,ps45i}(iPeak,iPlot));
                        tPKprog45 = abs(mean(dataToPlot1st{k,ps45i}(iPeak,iPlot))-mean(dataToPlot2nd{k,ps45i}(iPeak,iPlot)));
                        PKss = mean(dataToPlot1st{k,ssi}(iPeak,iPlot))-mean(dataToPlot2nd{k,ssi}(iPeak,iPlot));
                        tPKss = abs(mean(dataToPlot1st{k,ssi}(iPeak,iPlot))-mean(dataToPlot2nd{k,ssi}(iPeak,iPlot)));
                        trIndex{n,m}=(PKprog45-PKss)./(tPKprog45+tPKss);
                end
        end
        h(n,m) = notBoxPlot(trIndex{n,m},x(n,m));hold on
        set(h(n,m).data,'MarkerFaceColor',cReg(k,:))
        % save the maximum y values for significance plot alignment
        maxbl = nanmax(trIndex{n,m});
        if isempty(maxbl)
            maxbl = 0;
        end
        my(n) = max(my(n),maxbl);
        trIndex{n,m}(isnan(trIndex{n,m})) = [];
    end
    % Perform Mann-Whitney-Wilcoxon non parametric test and plot
    % significance indicators if p<0.05
    nck = nchoosek(1:m,2);
    %This is all done just to help with the plots
    nck = [nck; nck(end,:)];
    nph = size(nck,1) - 1; % number of comparisons made (used for adjusting the
    % minimum P value required for significance - a.k.a. Bonferroni correction)
    sameComparison = 1;
    for ii = 1 : size(nck,1) - 1
        if nck(ii+1) == nck(ii)
            sameComparison = sameComparison + 1;
        else
            sameComparison = 1;
        end
    end
    for ii = 1 : size(nck,1) - 1
        stats{ii} = mwwtest(trIndex{n,nck(ii,1)},trIndex{n,nck(ii,2)},0);
        % disp([compNames{n},'-',regions{neuronsToPlot(nck(ii,1))},'-',...
        %     regions{neuronsToPlot(nck(ii,2))},'-',num2str(stats{ii}.p(2))])
        if stats{ii}.p < 0.05/nph
            % uncomment below if you want to display the p-values per pair
            disp([compNames{n},'-',regions{neuronsToPlot(nck(ii,1))},'-',...
                regions{neuronsToPlot(nck(ii,2))},'-',num2str(stats{ii}.p(2))])
            sigline([x(n,nck(ii,1)),x(n,nck(ii,2))],[],my(n)+my(n)*0.05*ii);
        end
    end
end
limsX = xlim;
line(limsX,[0,0],'LineStyle','--','Color','k','LineWidth',1)
set(gca,'XTick',x(:,end)','XTickLabel',compNames,'FontSize',18,...
    'YLim',limsY)
set(gca,'FontSize',14,'YLim',limsY)
ylabel('Discrimination Index (peak)')

for m=1:length(neuronsToPlot)
    hold on
    l(m)=plot(nan,nan,'Color',cReg(neuronsToPlot(m),:),'LineWidth',4);
end
set(gca,'XLim',limsX)
lh=legend(l,regions(neuronsToPlot),'Box','off','FontSize',24,'Location',...
    'Best');
clearvars l lh h

%% Extended Data Figure 5a, 5c, 5d
% --------------------------- User Input ----------------------------------
% which OF stimuli to plot. Defined by their index in the Stimuli cell
% array. Plotted in the same order. Leave blank to plot all stimuli.
stimToPlot = [2 10 15];

% which neurons to plot. Defined by their index in the regions cell array.
% Plotted in the same order. Leave blank to plot all neurons.
neuronsToPlot = [8,4,11,5,12];

% Which calcium indicator(s) to plot. Defined by the genotype field in OF
% struct array. Leave blank to plot everything. Note that it is possible to
% pool data from multiple indicators (e.g. gcampsToPlot = {GC6f','GC7f'})
% Accepted arguments: {'GC6f','GC7f','sytGC7f','GC8f'}
% Default: gcampsToPlot = {};
gcampsToPlot = {'GC7f','sytGC7f'};

% Which prep types to plot. Defined by the prepType field in OF struct
% array. Leave blank to plot all experiments. Keep in mind that all of the
% OF sensitivity experiments in walking prep were done when the fly was
% stationary and not moving/grooming/struggling etc.
% Accepted arguments: ['walking','immobilized'] Default: prepToPlot = '';
prepToPlot = 'immobilized';

% Neurons from which hemisphere to plot. Defined by the somaSide field in
% OF struct array. Type 'Both' to plot all neurons. In that case, responses
% from Right hemisphere will be flipped to match the Left hemisphere neural
% responses and plotted together.
% Accepted arguments: ['L','R','Both'] Default: sideToPlot = 'Both';
sideToPlot = 'Both';

% specify normalization type. Type 'DFF' for (F-Fo)/Fo normalization.
% 'Max' will normalize the data per ROI to the maximum DFF response
% recorded in that ROI (see variable maxes).
% 'Zscore' will z-score the data per ROI based on the baseline mean and
% standard deviation for each trial done on that ROI.
% 'Znorm' normalizes Z-score data by the max Z-score recorded from each ROI
% Accepted arguments: ['DFF','Max','Zscore','Znorm']
normType = 'Zscore';

% Specify plot type. Type 'Line' for plotting each fly with a different
% colored line together with the mean trace. 'Shaded' will plot the std as
% a shaded area around the mean trace instead.
% Accepted arguments: ['Line','Shaded']
plotType = 'Shaded';

% Define Y-axis limits for each normalization type.
switch normType
    case 'Max'
        limsY = [-0.5 1.1]; % Y-axis limits
    case 'DFF'
        limsY = [-3 10]; % Y-axis limits
    case 'Zscore'
        switch plotType
            case 'Line'
                limsY = [-6 50]; % Y-axis limits
            case 'Shaded'
                limsY = [-3.5 28]; % Y-axis limits
        end
    case 'Znorm'
        limsY = [-0.4 1.05]; % Y-axis limits
end
%----------------------------end of user input-----------------------------
H2flip = 1;
H2i = find(contains(regions,'H2axon')); % index of H2
% Define the subset of data that will be plotted. This will be used in the
% plotSelectDFF function to plot only the specified subset of data as well
% as to calculate and display number of flies & ROIs being plotted
if isempty(stimToPlot)
    stimToPlot = 1:length(Stimuli);
end
if isempty(neuronsToPlot)
    neuronsToPlot = 1:length(regions);
end
prepIndices = cell(length(regions),1);
gcIndices = cell(length(regions),1);
sideIndices = cell(length(regions),1);
prepIndicesFly = cell(length(regions),1);
gcIndicesFly = cell(length(regions),1);
sideIndicesFly = cell(length(regions),1);
if ~isempty(prepToPlot)
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        if isempty(prepIndices{r})
            prepIndices{r} = ...
                contains(PrepTypePerROI{r},prepToPlot);
            prepIndicesFly{r} = ...
                contains(PrepType{r},prepToPlot);
        else
            prepIndices{r} = (prepIndices{r} | ...
                contains(PrepTypePerROI{r},prepToPlot));
            prepIndicesFly{r} = (prepIndicesFly{r} | ...
                contains(PrepType{r},prepToPlot));
        end
    end
else
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        prepIndices{r} = true(size(PrepTypePerROI{r}));
        prepIndicesFly{r} = true(size(PrepType{r}));
    end
end
if ~isempty(gcampsToPlot)
    for p = 1:length(gcampsToPlot) %for each selected GCaMP type
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            if isempty(gcIndices{r})
                gcIndices{r} = ...
                    ismember(GCaMPTypePerROI{r},gcampsToPlot{p});
                gcIndicesFly{r} = ...
                    ismember(GCaMPType{r},gcampsToPlot{p});
            else
                gcIndices{r} = (gcIndices{r} | ...
                    ismember(GCaMPTypePerROI{r},gcampsToPlot{p}));
                gcIndicesFly{r} = (gcIndicesFly{r} | ...
                    ismember(GCaMPType{r},gcampsToPlot{p}));
            end
        end
    end
else
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        gcIndices{r} = true(size(GCaMPTypePerROI{r}));
        gcIndicesFly{r} = true(size(GCaMPType{r}));
    end
end
switch sideToPlot
    case 'Both'
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            sideIndices{r} = true(size(SomaSidePerROI{r}));
            sideIndicesFly{r} = true(size(SomaSide{r}));
        end
        switch normType
            case 'Zscore'
                dataToPlot1st = OFRes_1st_Z_R;
                dataToPlot2nd = OFRes_2nd_Z_R;
                if H2flip
                    dataToPlot1st(H2i,:) = OFRes_1st_Z(H2i,:);
                    dataToPlot2nd(H2i,:) = OFRes_2nd_Z(H2i,:);
                end
            case 'Max'
                dataToPlot1st = OFRes_1st_normM;
                dataToPlot2nd = OFRes_2nd_normM;
            case 'DFF'
                dataToPlot1st = OFRes_1st;
                dataToPlot2nd = OFRes_2nd;
            case 'Znorm'
                dataToPlot1st = OFRes_1st_Znorm;
                dataToPlot2nd = OFRes_2nd_Znorm;
            otherwise
                error(['Normalization type not specified correctly. ',...
                    'Use normType=[''Max'',''Zscore'',''Znorm'' ',...
                    'or ''Raw'']'])
        end
    case {'L','R'}
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            if isempty(sideIndices{r})
                sideIndices{r} = ...
                    contains(SomaSidePerROI{r},sideToPlot);
                sideIndicesFly{r} = ...
                    contains(SomaSide{r},sideToPlot);
            else
                sideIndices{r} = (sideIndices{r} | ...
                    contains(SomaSidePerROI{r},sideToPlot));
                sideIndicesFly{r} = (sideIndicesFly{r} | ...
                    contains(SomaSide{r},sideToPlot));
            end
        end
        switch normType
            case 'Zscore'
                dataToPlot1st = OFRes_Raw_1st_Z;
                dataToPlot2nd = OFRes_Raw_2nd_Z;
            case 'Max'
                dataToPlot1st = OFRes_Raw_1st_normM;
                dataToPlot2nd = OFRes_Raw_2nd_normM;
            case 'Raw'
                dataToPlot1st = OFRes_Raw_1st;
                dataToPlot2nd = OFRes_Raw_2nd;
            case 'Znorm'
                dataToPlot1st = OFRes_1st_Znorm;
                dataToPlot2nd = OFRes_2nd_Znorm;
            otherwise
                error(['Normalization type not specified correctly. ',...
                    'Use normType=[''Max'',''Zscore'',''Znorm'' ',...
                    'or ''Raw'']'])
        end
    otherwise
        error(['Plotting side not specified correctly. ',...
            'Use sideToPlot = [''L'',''R'' or ''Both'']'])
end
%--------------------------------------------------------------------------
% Extended Data Figure 5a
f=plotSelectDFF(dataToPlot1st,dataToPlot2nd,neuronsToPlot,...
    stimToPlot,minFrameRate,cReg,prepIndices,prepIndicesFly,...
    gcIndices,gcIndicesFly,sideIndices,sideIndicesFly,Stimuli,...
    regions,limsY,plotType,normType);
%--------------------------------------------------------------------------
% Extended Data Figure 5c
stimToPlot = [17 16];
f=plotSelectDFF(dataToPlot1st,dataToPlot2nd,neuronsToPlot,...
    stimToPlot,minFrameRate,cReg,prepIndices,prepIndicesFly,...
    gcIndices,gcIndicesFly,sideIndices,sideIndicesFly,Stimuli,...
    regions,limsY,plotType,normType);
%--------------------------------------------------------------------------
% Extended Data Figure 5d
stimToPlot = [14 7];
f=plotSelectDFF(dataToPlot1st,dataToPlot2nd,neuronsToPlot,...
    stimToPlot,minFrameRate,cReg,prepIndices,prepIndicesFly,...
    gcIndices,gcIndicesFly,sideIndices,sideIndicesFly,Stimuli,...
    regions,limsY,plotType,normType);
%--------------------------------------------------------------------------
%% Extended Data Figure 5b
% --------------------------- User Input ----------------------------------
% which OF stimuli to plot. Defined by their index in the Stimuli cell
% array. Plotted in the same order. Leave blank to plot all stimuli.
stimToPlot = find(contains(StimuliSTPooled,'DotYaw_4pxsq_CCW'));

% which neurons to plot. Defined by their index in the regions cell array.
% Plotted in the same order. Leave blank to plot all neurons.
neuronsToPlot = [8,4,11,5,12];

% Which calcium indicator(s) to plot. Defined by the genotype field in OF
% struct array. Leave blank to plot everything. Note that it is possible to
% pool data from multiple indicators (e.g. gcampsToPlot = {GC6f','GC7f'})
% Accepted arguments: {'GC6f','GC7f','sytGC7f'} Default: gcampsToPlot = {};
gcampsToPlot = {'GC7f','sytGC7f'};

% Which prep types to plot. Defined by the prepType field in OF struct
% array. Leave blank to plot all experiments.
% Accepted arguments: ['walking','immobilized'] Default: prepToPlot = '';
prepToPlot = 'immobilized';

% Neurons from which hemisphere to plot. Defined by the somaSide field in
% OF struct array. Type 'Both' to plot all neurons. In that case, responses
% from Right hemisphere will be flipped to match the Left hemisphere neural
% responses and plotted together.
% Accepted arguments: ['L','R','Both'] Default: sideToPlot = 'Both';
sideToPlot = 'Both';

% specify normalization type. Type 'Raw' for no normalization.
% 'Max' will normalize the data per ROI to the maximum response recorded in
% that ROI (see variable maxes).
% 'Zscore' will zscore the data per ROI based on the mean and standard
% deviation of the entire recording session done on that ROI (see variables
% mu and sigma). 'Znorm' normalizes Z-score data by the max Z-score
% observed per fly per neuron
% Accepted arguments: ['Max','Zscore','Raw','Znorm'] Default: 'Raw';
normType = 'Zscore';

% Specify plot type. Type 'Line' for plotting each fly with a different
% colored line together with the mean trace. 'Shaded' will plot the std as
% a shaded area around the mean trace instead.
% Accepted arguments: ['Line','Shaded']
plotType = 'Shaded';

switch normType
    case 'Max'
        limsY = [-0.5 1.1]; % Y-axis limits
    case 'Raw'
        limsY = [-3 10]; % Y-axis limits
    case 'Zscore'
        switch plotType
            case 'Line'
                limsY = [-5 30]; % Y-axis limits
            case 'Shaded'
                %                 limsY = [-3.5 30]; % Y-axis limits
                %                 limsY = [-5 15]; % Y-axis limits
                limsY = [-5 35]; % Y-axis limits
        end
    case 'Znorm'
        limsY = [-0.25 1]; % Y-axis limits
end
%----------------------------end of user input-----------------------------

% Define the subset of data that will be plotted. This will be used in the
% plotSelectDFF function to plot only the specified subset of data as well
% as to calculate and display number of ROIs being plotted
if isempty(neuronsToPlot)
    neuronsToPlot = 1:length(regions);
end
prepIndices = cell(length(regions),1);
gcIndices = cell(length(regions),1);
sideIndices = cell(length(regions),1);
prepIndicesFly = cell(length(regions),1);
gcIndicesFly = cell(length(regions),1);
sideIndicesFly = cell(length(regions),1);
if ~isempty(prepToPlot)
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        if isempty(prepIndices{r})
            prepIndices{r} = ...
                contains(PrepTypePerROI{r},prepToPlot);
            prepIndicesFly{r} = ...
                contains(PrepType{r},prepToPlot);
        else
            prepIndices{r} = (prepIndices{r} | ...
                contains(PrepTypePerROI{r},prepToPlot));
            prepIndicesFly{r} = (prepIndicesFly{r} | ...
                contains(PrepType{r},prepToPlot));
        end
    end
else
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        prepIndices{r} = true(size(PrepTypePerROI{r}));
        prepIndicesFly{r} = true(size(PrepType{r}));
    end
end
if ~isempty(gcampsToPlot)
    for p = 1:length(gcampsToPlot) %for each selected GCaMP type
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            if isempty(gcIndices{r})
                gcIndices{r} = ...
                    ismember(GCaMPTypePerROI{r},gcampsToPlot{p});
                gcIndicesFly{r} = ...
                    ismember(GCaMPType{r},gcampsToPlot{p});
            else
                gcIndices{r} = (gcIndices{r} | ...
                    ismember(GCaMPTypePerROI{r},gcampsToPlot{p}));
                gcIndicesFly{r} = (gcIndicesFly{r} | ...
                    ismember(GCaMPType{r},gcampsToPlot{p}));
            end
        end
    end
else
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        gcIndices{r} = true(size(GCaMPTypePerROI{r}));
        gcIndicesFly{r} = true(size(GCaMPType{r}));
    end
end
switch sideToPlot
    case 'Both'
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            sideIndices{r} = true(size(SomaSidePerROI{r}));
            sideIndicesFly{r} = true(size(SomaSide{r}));
        end
        switch normType
            case 'Zscore'
                dataToPlot = STResP_Z;
            case 'Max'
                dataToPlot = STResP_normM;
            case 'Raw'
                dataToPlot = STResP;
            case 'Znorm'
                %                 dataToPlot1st = OFRes_1st_Znorm;
                %                 dataToPlot2nd = OFRes_2nd_Znorm;
            otherwise
                error(['Normalization type not specified correctly. ',...
                    'Use normType=[''Max'',''Zscore'',''Znorm'' ',...
                    'or ''Raw'']'])
        end
    case {'L','R'}
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            if isempty(sideIndices{r})
                sideIndices{r} = ...
                    contains(SomaSidePerROI{r},sideToPlot);
                sideIndicesFly{r} = ...
                    contains(SomaSide{r},sideToPlot);
            else
                sideIndices{r} = (sideIndices{r} | ...
                    contains(SomaSidePerROI{r},sideToPlot));
                sideIndicesFly{r} = (sideIndicesFly{r} | ...
                    contains(SomaSide{r},sideToPlot));
            end
        end
        switch normType
            case 'Zscore'
                dataToPlot = STResP_Raw_Z;
            case 'Max'
                dataToPlot = STResP_Raw_normM;
            case 'Raw'
                dataToPlot1st = STResP_Raw;
            case 'Znorm'
                %                 dataToPlot1st = OFRes_1st_Znorm;
                %                 dataToPlot2nd = OFRes_2nd_Znorm;
            otherwise
                error(['Normalization type not specified correctly. ',...
                    'Use normType=[''Max'',''Zscore'',''Znorm'' ',...
                    'or ''Raw'']'])
        end
    otherwise
        error(['Plotting side not specified correctly. ',...
            'Use sideToPlot = [''L'',''R'' or ''Both'']'])
end
%--------------------------------------------------------------------------
% Extended Data Figure 5b z-scored traces per ST stimulus (Left)
f=plotSelectDFFST(dataToPlot,neuronsToPlot,...
    stimToPlot,minFrameRate,cReg,prepIndices,prepIndicesFly,...
    gcIndices,gcIndicesFly,sideIndices,sideIndicesFly,StimuliSTPooled,...
    regions,limsY,plotType);

%--------------------------------------------------------------------------
% Extended Data Figure 5b speed tuning curves (Right)
neuronsToPlot = [8,4,11,5,12];
PlotType = 'bar'; % accepted arguments  = {'bar','shaded'}

% initialize speed tuning curve related variables
STcurves = cell(length(regions),1);
STcurves_Z = cell(length(regions),1);

% create a new figure
f=figure('WindowState','maximized');
hold on
for i =  neuronsToPlot % for each unique neuron-region combination
    switch regions{i} % define which OF ST curve will be plotted per neuron
        case {'H2axon'}
            stimIndex = find(contains(StimuliSTPooled,'DotYaw_4pxsq_CW'));
        otherwise
            stimIndex = find(contains(StimuliSTPooled,'DotYaw_4pxsq_CCW'));
    end
    stimNames = StimuliSTPooled(stimIndex);

    % initialize ST curve related variables
    speeds = cell(1,length(stimNames));
    STcurves{i} = nan(size(STResP_Z{i,stimIndex(1)},2),length(stimNames));

    % Define the indices for each trial. Used to define the range for auc
    % calculations (seconds 1-2 for baseline, 2-6 for stimulus-evoked auc)
    frameNum = size(STResP_Z{i,1},1);
    TimeIndices = linspace(0,frameNum/(minFrameRate/1000),frameNum);

    for spd=1:length(stimNames) % for each stimulus speed
        % get the value of speeds used for the x-axis label
        dashes = strfind(stimNames{spd},'_');
        dashes2 = strfind(stimNames{spd},'px/s');
        tempspeed = stimNames{spd}(dashes(end)+1:dashes2(end)-1);
        speeds{spd} = str2double(tempspeed)*2.25;
        if ~isempty(STResP_Z{i,stimIndex(spd)}) % if there is imaging data
            if ~isempty(STResP_Z{i,stimIndex(spd)}...
                    (~isnan(STResP_Z{i,stimIndex(spd)})))%if data has value
                % auc of the stimulus evoked calcium response (seconds 2-6)
                STcurves{i}(:,spd) = trapz(STResP_Z{i,stimIndex(spd)}...
                    ((find(TimeIndices<2,1,'last')+1):end,:));
                STcurves{i}(:,spd) = trapz(STResP_Z{i,stimIndex(spd)}...
                    ((find(TimeIndices<2,1,'last')+1):...
                    (find(TimeIndices<4,1,'last')),:));
            end
        end
    end
    STcurves{i}(isnan(STcurves{i}(:,spd)),:) = []; % remove NaN data points
    % normalize each ROI based on its absolute maximum auc response
    STcurves_Z{i} = STcurves{i}./max(abs(STcurves{i}),[],2);
    switch PlotType
        case 'shaded'
            % plot the curves with shaded error bars
            shadedErrorBar(1:length(speeds),mean(STcurves_Z{i}),...
                std(STcurves_Z{i})/sqrt(size(STcurves_Z{i},1)),'lineprops',...
                {'Color',cReg(i,:),'LineWidth',3},'transparent',1)
        case 'bar'
            % plot the curves with line error bars
            eh = errorbar(1:length(speeds),mean(STcurves_Z{i}),...
                std(STcurves_Z{i})/sqrt(size(STcurves_Z{i},1)));
            set(eh, 'Color', cReg(i,:),'LineWidth',3, 'CapSize', 15)
        otherwise
            error('unrecognized PlotType variable')
    end
end

limsX = [0 length(speeds)+1];
% modify X-axis labels to match the speeds presented, add y-label
set(gca,'XTick',1:length(speeds),'XTickLabel',speeds,'XLim',...
    limsX,'FontSize',18)
xlabel('Stimulus speed (deg/s)');
ylabel('Normalized response');
% add a dashed line at zero for reference
line(limsX,[0,0],'LineStyle','--','Color','k','LineWidth',1)
% add imaginary lines for each neuron (used to make the legend)
l = gobjects(1,length(neuronsToPlot));
for k=1:length(neuronsToPlot)
    hold on
    l(k)=plot(nan,nan,'Color',cReg(neuronsToPlot(k),:),'LineWidth',4);
end
lh=legend(l,regions(neuronsToPlot),'Box','off','FontSize',24); % add legend
clearvars l lh eh
%--------------------------------------------------------------------------
%% Extended Data Figure 6b-f
% --------------------------- User Input ----------------------------------
neuronsToPlot = [8,4,11,5,12];
% Specify plot type. Type 'Line' for plotting each fly with a different
% colored line together with the mean trace. 'Shaded' will plot the std as
% a shaded area around the mean trace instead.
% Accepted arguments: ['Line','Shaded']
plotType = 'Shaded';
%----------------------------end of user input-----------------------------

% Translational stimuli used (as of 210421)
StimTR = {'Sideslip','ProgSide315','Progressive','ProgSide45'};

% Find the corresponding indices for the stimuli defined in StimTR in the
% Stimuli cell array
iTR = nan(1,length(StimTR));
for i = 1:length(iTR)
    iTR(i)=find(contains(Stimuli,StimTR(i))&~contains(Stimuli,'Left') & ...
        ~contains(Stimuli,'Right'));
end

limsY = [-3.5 30]; % Y-axis limits for HS
limsY2 = [-3.5 13]; % Y-axis limits for the rest
cRect = [0 0.5 0 0.1]; % color of the rectangle depicting stimulus motion

for i = 1:length(neuronsToPlot) % for each unique neuron-region combination
    f=figure('WindowState','maximized'); hold on
    n = neuronsToPlot(i);
    if contains(regions{n},'HS')
        lmY = limsY;
    else
        lmY = limsY2;
    end
    frameNum = size(OFRes_1st_Z_R{n,1},1);
    TimeAxisPD = linspace(0,frameNum/(minFrameRate/1000),frameNum);
    %     limsX = [min(TimeAxisPD) ceil(max(TimeAxisPD))];
    limsX = [min(TimeAxisPD)+1 ceil(max(TimeAxisPD))-0.5];

    for j = 1:length(StimTR) % for each OF stimulus
        if ~isempty(OFRes_1st_Z_R{n,iTR(j)}) % if there is imaging data
            CounterTRn = sum(~isnan(OFRes_1st_Z_R{n,iTR(j)}(1,:)));
            if contains(regions{n},'H2axon') % Flip to plot Left H2
                switch StimTR{j}
                    case 'Sideslip'
                        subplot(3,3,4)
                    case 'ProgSide315'
                        subplot(3,3,3)
                    case 'Progressive'
                        subplot(3,3,2)
                    case 'ProgSide45'
                        subplot(3,3,1)
                end
            else
                switch StimTR{j}
                    case 'Sideslip'
                        subplot(3,3,6)
                    case 'ProgSide315'
                        subplot(3,3,1)
                    case 'Progressive'
                        subplot(3,3,2)
                    case 'ProgSide45'
                        subplot(3,3,3)
                end
            end
            % plot individual neurites, and mean across neurites, 1st epoch
            switch plotType
                case 'Line'
                    plot(TimeAxisPD,OFRes_1st_Z_R{n,iTR(j)}); hold on;
                    plot(TimeAxisPD,nanmean(OFRes_1st_Z_R{n,iTR(j)},2),...
                        'Color',cReg(n,:),'LineWidth',4);
                case 'Shaded'
                    h=shadedErrorBar(TimeAxisPD,...
                        nanmean(OFRes_1st_Z_R{n,iTR(j)},2),...
                        nanstd(OFRes_1st_Z_R{n,iTR(j)},0,2)/...
                        sqrt(CounterTRn));hold on
                    set(h.mainLine,'Color',cReg(n,:),...
                        'LineWidth',4)
                    set(h.patch,'FaceColor',cReg(n,:),...
                        'FaceAlpha',0.2)
                    set(h.edge(1),'Color',cReg(n,:))
                    set(h.edge(2),'Color',cReg(n,:))
                    clearvars h
            end

            rectangle('Position',[2 lmY(1) 2 lmY(2)-lmY(1)],...
                'FaceColor', cRect);
            set(gca,'YLim',lmY,'XLim',limsX)
            switch StimTR{j}
                case 'ProgSide315'
                    ylabel('Response (z-score)')
                    xlabel('Time (sec)')
                    set(gca,'FontSize',18)
                otherwise
                    axis off
            end
            ltext = ['n=' num2str(CounterTRn)];
            th = text(0,0.9*lmY(2),ltext);
            set(th,'FontSize',18)
            line(limsX,[0 0],'Color','k','LineWidth',0.5,'LineStyle','--')
        end

        if  ~isempty(OFRes_2nd_Z_R{n,iTR(j)}) % if there is imaging data
            if contains(regions{n},'H2axon') % Flip to plot Left H2
                switch StimTR{j}
                    case 'Sideslip'
                        subplot(3,3,6)
                    case 'ProgSide315'
                        subplot(3,3,7)
                    case 'Progressive'
                        subplot(3,3,8)
                    case 'ProgSide45'
                        subplot(3,3,9)
                end
            else
                switch StimTR{j}
                    case 'Sideslip'
                        subplot(3,3,4)
                    case 'ProgSide315'
                        subplot(3,3,9)
                    case 'Progressive'
                        subplot(3,3,8)
                    case 'ProgSide45'
                        subplot(3,3,7)
                end
            end
            % plot individual neurites, and mean across neurites, 2nd epoch
            switch plotType
                case 'Line'
                    plot(TimeAxisPD,OFRes_2nd_Z_R{n,iTR(j)}); hold on;
                    plot(TimeAxisPD,nanmean(OFRes_2nd_Z_R{n,iTR(j)},2),...
                        'color',cReg(n,:),'LineWidth',4);
                case 'Shaded'
                    h=shadedErrorBar(TimeAxisPD,...
                        nanmean(OFRes_2nd_Z_R{n,iTR(j)},2),...
                        nanstd(OFRes_2nd_Z_R{n,iTR(j)},0,2)/...
                        sqrt(CounterTRn));hold on
                    set(h.mainLine,'Color',cReg(n,:),...
                        'LineWidth',4)
                    set(h.patch,'FaceColor',cReg(n,:),...
                        'FaceAlpha',0.2)
                    set(h.edge(1),'Color',cReg(n,:))
                    set(h.edge(2),'Color',cReg(n,:))
                    clearvars h
            end
            rectangle('Position',[2 lmY(1) 2 lmY(2)-lmY(1)],...
                'FaceColor', cRect);
            set(gca,'YLim',lmY,'XLim',limsX)
            axis off
            ltext = ['n=' num2str(CounterTRn)];
            th = text(0,0.9*lmY(2),ltext);
            set(th,'FontSize',18)
            line(limsX,[0 0],'Color','k','LineWidth',0.5,'LineStyle','--')
        end
    end
    % Add number of neurons to the middle subplot
    subplot(3,3,5)
    ltext = [regions{n}];
    th = text(0.3,0.5,ltext);
    set(th,'FontSize',22)
    box off
    axis off
end

%% Extended Data Figure 6g-h

% --------------------------- User Input ----------------------------------
% which neurons to plot. Defined by their index in the regions cell array.
% Plotted in the same order. Leave blank to plot all neurons.

neuronsToPlot = [11,8,4,5,12];

% Which calcium indicator(s) to plot. Defined by the genotype field in OF
% struct array. Leave blank to plot everything. Note that it is possible to
% pool data from multiple indicators (e.g. gcampsToPlot = {GC6f','GC7f'})
% Accepted arguments: {'GC6f','GC7f','sytGC7f'} Default: gcampsToPlot = {};
gcampsToPlot = {'sytGC7f','GC7f'};

% Which prep types to plot. Defined by the prepType field in OF struct
% array. Leave blank to plot all experiments.
% Accepted arguments: ['walking','immobilized'] Default: prepToPlot = '';
prepToPlot = 'immobilized';

% Neurons from which hemisphere to plot. Defined by the somaSide field in
% OF struct array. Type 'Both' to plot all neurons. In that case, responses
% from Right hemisphere will be flipped to match the Left hemisphere neural
% responses and plotted together.
% Accepted arguments: ['L','R','Both'] Default: sideToPlot = 'Both';
sideToPlot = 'Both';

% specify normalization type. Type 'Raw' for no normalization.
% 'Max' will normalize the data per ROI to the maximum response recorded in
% that ROI (see variable maxes).
% 'Zscore' will zscore the data per ROI based on the mean and standard
% deviation of the entire recording session done on that ROI (see variables
% mu and sigma). 'Znorm' normalizes Z-score data by the max Z-score
% observed per fly per neuron
% Accepted arguments: ['Zscore','Raw'] Default: 'Raw';
normType = 'Zscore';

% Specify plot type. Type 'Box' for plotting each fly with a different
% colored dot together with the population mean and sd as boxes.
% Currently only 'Box' is supported.
% Accepted arguments: ['Box']
plotType = 'Box';

% Specify grouping. Type 'neuron' for plotting all stimuli for each neuron
% grouped together. 'stimulus' will group all neurons for each stimulus
% Accepted arguments: ['neuron','stimulus']
groupType = 'neuron';

%----------------------------end of user input-----------------------------

peakbin = 500; % time bin for the averaged peak response (ms)
H2Flip = 0;
iH2 = find(contains(regions,'H2axon')); % index of H2
% index for recordings that will come from the right hemisphere
[~,loc] = ismember(StimuliFlipped,Stimuli);

clearvars l lh h stats
stm = {'Yaw_','Progressive','YawRightFlipped','RollInv','Pitch'};
auci = find(contains(Stimuli,stm) & BinocularIndices);
yi = auci(contains(Stimuli(auci),'Yaw_'));
yfi = auci(contains(Stimuli(auci),'YawRightFlipped'));
pri = auci(contains(Stimuli(auci),'Progressive'));
ssi = find(contains(Stimuli,'Sideslip'));
ps45i = find(contains(Stimuli,'ProgSide45'));

% Define the subset of data that will be plotted. This will be used in the
% plotSelectAUC function to plot only the specified subset of data as well
% as to calculate and display number of flies & ROIs being plotted

if isempty(neuronsToPlot)
    neuronsToPlot = 1:length(regions);
end

prepIndices = cell(length(regions),1);
gcIndices = cell(length(regions),1);
sideIndices = cell(length(regions),1);
prepIndicesFly = cell(length(regions),1);
gcIndicesFly = cell(length(regions),1);
sideIndicesFly = cell(length(regions),1);

if ~isempty(prepToPlot)
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        if isempty(prepIndices{r})
            prepIndices{r} = ...
                contains(PrepTypePerROI{r},prepToPlot);
            prepIndicesFly{r} = ...
                contains(PrepType{r},prepToPlot);
        else
            prepIndices{r} = (prepIndices{r} | ...
                contains(PrepTypePerROI{r},prepToPlot));
            prepIndicesFly{r} = (prepIndicesFly{r} | ...
                contains(PrepType{r},prepToPlot));
        end
    end
else
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        prepIndices{r} = true(size(PrepTypePerROI{r}));
        prepIndicesFly{r} = true(size(PrepType{r}));
    end
end
if ~isempty(gcampsToPlot)
    for p = 1:length(gcampsToPlot) %for each selected GCaMP type
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            if isempty(gcIndices{r})
                gcIndices{r} = ...
                    ismember(GCaMPTypePerROI{r},gcampsToPlot{p});
                gcIndicesFly{r} = ...
                    ismember(GCaMPType{r},gcampsToPlot{p});
            else
                gcIndices{r} = (gcIndices{r} | ...
                    ismember(GCaMPTypePerROI{r},gcampsToPlot{p}));
                gcIndicesFly{r} = (gcIndicesFly{r} | ...
                    ismember(GCaMPType{r},gcampsToPlot{p}));
            end
        end
    end
else
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        gcIndices{r} = true(size(GCaMPTypePerROI{r}));
        gcIndicesFly{r} = true(size(GCaMPType{r}));
    end
end
switch sideToPlot
    case 'Both'
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            sideIndices{r} = true(size(SomaSidePerROI{r}));
            sideIndicesFly{r} = true(size(SomaSide{r}));
        end
        switch normType
            case 'Zscore'
                dataToPlot1st = OFRes_1st_Z;
                dataToPlot2nd = OFRes_2nd_Z;
            case 'Max'
                dataToPlot1st = Auc_1st_normM;
                dataToPlot2nd = Auc_2nd_normM;
            case 'Raw'
                dataToPlot1st = Auc_1st_bs;
                dataToPlot2nd = Auc_2nd_bs;
            otherwise
                error(['Normalization type not specified correctly. ',...
                    'Use normType=[''Max'',''Zscore'', ',...
                    'or ''Raw'']'])
        end
    case {'L','R'}
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            if isempty(sideIndices{r})
                sideIndices{r} = ...
                    contains(SomaSidePerROI{r},sideToPlot);
                sideIndicesFly{r} = ...
                    contains(SomaSide{r},sideToPlot);
            else
                sideIndices{r} = (sideIndices{r} | ...
                    contains(SomaSidePerROI{r},sideToPlot));
                sideIndicesFly{r} = (sideIndicesFly{r} | ...
                    contains(SomaSide{r},sideToPlot));
            end
        end
        switch normType
            case 'Zscore'
                dataToPlot1st = Auc_1st_bs_Z;
                dataToPlot2nd = Auc_2nd_bs_Z;
            case 'Max'
                dataToPlot1st = Auc_1st_normM;
                dataToPlot2nd = Auc_2nd_normM;
            case 'Raw'
                dataToPlot1st = Auc_1st_bs;
                dataToPlot2nd = Auc_2nd_bs;
            otherwise
                error(['Normalization type not specified correctly. ',...
                    'Use normType=[''Max'',''Zscore'', ',...
                    'or ''Raw'']'])
        end
    otherwise
        error(['Plotting side not specified correctly. ',...
            'Use sideToPlot = [''L'',''R'' or ''Both'']'])
end

frameNum = size(dataToPlot1st{1,1},1);
TimeAxisPD = linspace(0,frameNum/(minFrameRate/1000),frameNum);
iPeak = find(TimeAxisPD>4-peakbin/1000,1,'first'):...
    find(TimeAxisPD<4,1,'last');

compNames = {'Prog-Sideslip','ProgSide45-Sideslip'};
ncomps = length(compNames);
sp=0.5;
x = zeros(ncomps,size(neuronsToPlot,2));
trIndex = cell(ncomps,size(neuronsToPlot,2));
f=figure('WindowState','maximized');
limsY = [-1.2 1.7];
my = zeros(1,ncomps);

for n = 1:ncomps
    for k = neuronsToPlot
        m = find(neuronsToPlot==k);
        x(n,m) = (n-1)*length(neuronsToPlot) + 1 + (m-1)*sp;
        iPlot = prepIndices{k} & gcIndices{k} & sideIndices{k};
        switch compNames{n}
            case 'Prog-Yaw'
                switch regions{k}
                    case {'H2axon'}
                        PKprog = mean(dataToPlot2nd{k,pri}(iPeak,iPlot))-mean(dataToPlot1st{k,pri}(iPeak,iPlot));
                        tPKprog = abs(mean(dataToPlot2nd{k,pri}(iPeak,iPlot))-mean(dataToPlot1st{k,pri}(iPeak,iPlot)));
                        PKyaw = mean(dataToPlot1st{k,yi}(iPeak,iPlot))-mean(dataToPlot2nd{k,yi}(iPeak,iPlot));
                        tPKyaw = abs(mean(dataToPlot1st{k,yi}(iPeak,iPlot))-mean(dataToPlot2nd{k,yi}(iPeak,iPlot)));
                        trIndex{n,m}=(PKprog-PKyaw)./(tPKprog+tPKyaw);
                    otherwise
                        PKprog = mean(dataToPlot1st{k,pri}(iPeak,iPlot))-mean(dataToPlot2nd{k,pri}(iPeak,iPlot));
                        tPKprog = abs(mean(dataToPlot1st{k,pri}(iPeak,iPlot))-mean(dataToPlot2nd{k,pri}(iPeak,iPlot)));
                        PKyaw = mean(dataToPlot2nd{k,yi}(iPeak,iPlot))-mean(dataToPlot1st{k,yi}(iPeak,iPlot));
                        tPKyaw = abs(mean(dataToPlot2nd{k,yi}(iPeak,iPlot))-mean(dataToPlot1st{k,yi}(iPeak,iPlot)));
                        trIndex{n,m}=(PKprog-PKyaw)./(tPKprog+tPKyaw);
                end
            case 'YawFlipped-Yaw'
                switch regions{k}
                    case {'H2axon'}
                        PKyf = mean(dataToPlot1st{k,yfi}(iPeak,iPlot))-mean(dataToPlot2nd{k,yfi}(iPeak,iPlot));
                        tPKyf = abs(mean(dataToPlot2nd{k,yfi}(iPeak,iPlot))-mean(dataToPlot1st{k,yfi}(iPeak,iPlot)));
                        PKyaw = mean(dataToPlot1st{k,yi}(iPeak,iPlot))-mean(dataToPlot2nd{k,yi}(iPeak,iPlot));
                        tPKyaw = abs(mean(dataToPlot1st{k,yi}(iPeak,iPlot))-mean(dataToPlot2nd{k,yi}(iPeak,iPlot)));
                        trIndex{n,m}=(PKyf-PKyaw)./(tPKyf+tPKyaw);
                    case {'H2rnaxon'}
                        PKyf = mean(dataToPlot1st{k,yfi}(iPeak,iPlot))-mean(dataToPlot2nd{k,yfi}(iPeak,iPlot));
                        tPKyf = abs(mean(dataToPlot2nd{k,yfi}(iPeak,iPlot))-mean(dataToPlot1st{k,yfi}(iPeak,iPlot)));
                        PKyaw = mean(dataToPlot2nd{k,yi}(iPeak,iPlot))-mean(dataToPlot1st{k,yi}(iPeak,iPlot));
                        tPKyaw = abs(mean(dataToPlot1st{k,yi}(iPeak,iPlot))-mean(dataToPlot2nd{k,yi}(iPeak,iPlot)));
                        trIndex{n,m}=(PKyf-PKyaw)./(tPKyf+tPKyaw);
                    otherwise
                        PKyf = mean(dataToPlot2nd{k,yfi}(iPeak,iPlot))-mean(dataToPlot1st{k,yfi}(iPeak,iPlot));
                        tPKyf = abs(mean(dataToPlot1st{k,yfi}(iPeak,iPlot))-mean(dataToPlot2nd{k,yfi}(iPeak,iPlot)));
                        PKyaw = mean(dataToPlot2nd{k,yi}(iPeak,iPlot))-mean(dataToPlot1st{k,yi}(iPeak,iPlot));
                        tPKyaw = abs(mean(dataToPlot2nd{k,yi}(iPeak,iPlot))-mean(dataToPlot1st{k,yi}(iPeak,iPlot)));
                        trIndex{n,m}=(PKyf-PKyaw)./(tPKyf+tPKyaw);
                end

            case 'Prog-Sideslip'
                switch regions{k}
                    case {'H2axon'}
                        PKss = mean(dataToPlot2nd{k,ssi}(iPeak,iPlot))-mean(dataToPlot1st{k,ssi}(iPeak,iPlot));
                        tPKss = abs(mean(dataToPlot2nd{k,ssi}(iPeak,iPlot))-mean(dataToPlot1st{k,ssi}(iPeak,iPlot)));
                        PKprog = mean(dataToPlot2nd{k,pri}(iPeak,iPlot))-mean(dataToPlot1st{k,pri}(iPeak,iPlot));
                        tPKprog = abs(mean(dataToPlot2nd{k,pri}(iPeak,iPlot))-mean(dataToPlot1st{k,pri}(iPeak,iPlot)));
                        trIndex{n,m}=(PKprog-PKss)./(tPKprog+tPKss);
                    otherwise
                        PKss = mean(dataToPlot1st{k,ssi}(iPeak,iPlot))-mean(dataToPlot2nd{k,ssi}(iPeak,iPlot));
                        tPKss = abs(mean(dataToPlot1st{k,ssi}(iPeak,iPlot))-mean(dataToPlot2nd{k,ssi}(iPeak,iPlot)));
                        PKprog = mean(dataToPlot1st{k,pri}(iPeak,iPlot))-mean(dataToPlot2nd{k,pri}(iPeak,iPlot));
                        tPKprog = abs(mean(dataToPlot1st{k,pri}(iPeak,iPlot))-mean(dataToPlot2nd{k,pri}(iPeak,iPlot)));
                        trIndex{n,m}=(PKprog-PKss)./(tPKprog+tPKss);
                end

            case 'ProgSide45-Sideslip'
                switch regions{k}
                    case {'H2axon'}
                        PKprog45 = mean(dataToPlot2nd{k,ps45i}(iPeak,iPlot))-mean(dataToPlot1st{k,ps45i}(iPeak,iPlot));
                        tPKprog45 = abs(mean(dataToPlot2nd{k,ps45i}(iPeak,iPlot))-mean(dataToPlot1st{k,ps45i}(iPeak,iPlot)));
                        PKss = mean(dataToPlot2nd{k,ssi}(iPeak,iPlot))-mean(dataToPlot1st{k,ssi}(iPeak,iPlot));
                        tPKss = abs(mean(dataToPlot2nd{k,ssi}(iPeak,iPlot))-mean(dataToPlot1st{k,ssi}(iPeak,iPlot)));
                        trIndex{n,m}=(PKprog45-PKss)./(tPKprog45+tPKss);
                    case {'H2rnaxon'}
                        PKprog45 = mean(dataToPlot1st{k,ps45i}(iPeak,iPlot))-mean(dataToPlot2nd{k,ps45i}(iPeak,iPlot));
                        tPKprog45 = abs(mean(dataToPlot1st{k,ps45i}(iPeak,iPlot))-mean(dataToPlot2nd{k,ps45i}(iPeak,iPlot)));
                        PKss = mean(dataToPlot1st{k,ssi}(iPeak,iPlot))-mean(dataToPlot2nd{k,ssi}(iPeak,iPlot));
                        tPKss = abs(mean(dataToPlot1st{k,ssi}(iPeak,iPlot))-mean(dataToPlot2nd{k,ssi}(iPeak,iPlot)));
                        trIndex{n,m}=(PKprog45-PKss)./(tPKprog45+tPKss);
                    otherwise
                        PKprog45 = mean(dataToPlot1st{k,ps45i}(iPeak,iPlot))-mean(dataToPlot2nd{k,ps45i}(iPeak,iPlot));
                        tPKprog45 = abs(mean(dataToPlot1st{k,ps45i}(iPeak,iPlot))-mean(dataToPlot2nd{k,ps45i}(iPeak,iPlot)));
                        PKss = mean(dataToPlot1st{k,ssi}(iPeak,iPlot))-mean(dataToPlot2nd{k,ssi}(iPeak,iPlot));
                        tPKss = abs(mean(dataToPlot1st{k,ssi}(iPeak,iPlot))-mean(dataToPlot2nd{k,ssi}(iPeak,iPlot)));
                        trIndex{n,m}=(PKprog45-PKss)./(tPKprog45+tPKss);
                end
        end
        h(n,m) = notBoxPlot(trIndex{n,m},x(n,m));hold on
        set(h(n,m).data,'MarkerFaceColor',cReg(k,:))
        % save the maximum y values for significance plot alignment
        maxbl = nanmax(trIndex{n,m});
        if isempty(maxbl)
            maxbl = 0;
        end
        my(n) = max(my(n),maxbl);
        trIndex{n,m}(isnan(trIndex{n,m})) = [];
    end
    % Perform Mann-Whitney-Wilcoxon non parametric test and plot
    % significance indicators if p<0.05
    nck = nchoosek(1:m,2);
    %This is all done just to help with the plots
    nck = [nck; nck(end,:)];
    nph = size(nck,1) - 1; % number of comparisons made (used for adjusting the
    % minimum P value required for significance - a.k.a. Bonferroni correction)
    sameComparison = 1;
    for ii = 1 : size(nck,1) - 1
        if nck(ii+1) == nck(ii)
            sameComparison = sameComparison + 1;
        else
            sameComparison = 1;
        end
    end
    for ii = 1 : size(nck,1) - 1
        stats{ii} = mwwtest(trIndex{n,nck(ii,1)},trIndex{n,nck(ii,2)},0);
        % disp([compNames{n},'-',regions{neuronsToPlot(nck(ii,1))},'-',...
        %     regions{neuronsToPlot(nck(ii,2))},'-',num2str(stats{ii}.p(2))])
        if stats{ii}.p < 0.05/nph
            % uncomment below if you want to display the p-values per pair
            disp([compNames{n},'-',regions{neuronsToPlot(nck(ii,1))},'-',...
                regions{neuronsToPlot(nck(ii,2))},'-',num2str(stats{ii}.p(2))])
            sigline([x(n,nck(ii,1)),x(n,nck(ii,2))],[],my(n)+my(n)*0.05*ii);
        end
    end
end
limsX = xlim;
line(limsX,[0,0],'LineStyle','--','Color','k','LineWidth',1)
set(gca,'XTick',x(:,end)','XTickLabel',compNames,'FontSize',18,...
    'YLim',limsY)
set(gca,'FontSize',14,'YLim',limsY)
ylabel('Discrimination Index (peak)')

for m=1:length(neuronsToPlot)
    hold on
    l(m)=plot(nan,nan,'Color',cReg(neuronsToPlot(m),:),'LineWidth',4);
end
set(gca,'XLim',limsX)
lh=legend(l,regions(neuronsToPlot),'Box','off','FontSize',24,'Location',...
    'Best');
clearvars l lh h

%% Extended Data Figure 8a-b

% Which compound OF stimuli to plot. Defined by their index in the Stimuli
% cell array. Plotted in the same order. Leave blank to plot all stimuli.

stimToPlot=find(contains(Stimuli,{'Flicker'})&contains(Stimuli,{'5px'}));

% which neurons to plot. Defined by their index in the regions cell array.
% Plotted in the same order. Leave blank to plot all neurons.

neuronsToPlot = [6,7,9,10,1,2]; % Figure 5b neurons

% Which calcium indicator(s) to plot. Defined by the genotype field in OF
% struct array. Leave blank to plot everything. Note that it is possible to
% pool data from multiple indicators (e.g. gcampsToPlot = {GC6f','GC7f'})
% Accepted arguments: {'GC6f','GC7f','sytGC7f'} Default: gcampsToPlot = {};
gcampsToPlot = {'GC6f'};

% Which prep types to plot. Defined by the prepType field in OF struct
% array. Leave blank to plot all experiments.
% Accepted arguments: ['walking','immobilized'] Default: prepToPlot = '';
prepToPlot = '';

% Neurons from which hemisphere to plot. Defined by the somaSide field in
% OF struct array. Type 'Both' to plot all neurons. In that case, responses
% from Right hemisphere will be flipped to match the Left hemisphere neural
% responses and plotted together.
% Accepted arguments: ['L','R','Both'] Default: sideToPlot = 'Both';
sideToPlot = 'Both';

% specify normalization type. Type 'Raw' for no normalization.
% 'Max' will normalize the data per ROI to the maximum response recorded in
% that ROI (see variable maxes).
% 'Zscore' will zscore the data per ROI based on the mean and standard
% deviation of the entire recording session done on that ROI (see variables
% mu and sigma). 'Znorm' normalizes Z-score data by the max Z-score
% observed per fly per neuron
% Accepted arguments: ['Max','Zscore','Raw','Znorm'] Default: 'Raw';
normType = 'Zscore';

% Specify plot type. Type 'Line' for plotting each fly with a different
% colored line together with the mean trace. 'Shaded' will plot the std as
% a shaded area around the mean trace instead.
% Accepted arguments: ['Line','Shaded']
plotType = 'Shaded';

switch normType
    case 'Max'
        limsY = [-0.5 1.1]; % Y-axis limits
    case 'Raw'
        limsY = [-3 10]; % Y-axis limits
    case 'Zscore'
        switch plotType
            case 'Line'
                limsY = [-5 30]; % Y-axis limits
            case 'Shaded'
                limsY = [-2.5 18]; % Y-axis limits
        end
    case 'Znorm'
        limsY = [-0.25 1]; % Y-axis limits
end
%----------------------------end of user input-----------------------------

% Define the subset of data that will be plotted. This will be used in the
% plotSelectDFF function to plot only the specified subset of data as well
% as to calculate and display number of flies & ROIs being plotted
if isempty(stimToPlot)
    stimToPlot = 1:length(Stimuli);
end
if isempty(neuronsToPlot)
    neuronsToPlot = 1:length(regions);
end

prepIndices = cell(length(regions),1);
gcIndices = cell(length(regions),1);
sideIndices = cell(length(regions),1);
prepIndicesFly = cell(length(regions),1);
gcIndicesFly = cell(length(regions),1);
sideIndicesFly = cell(length(regions),1);

if ~isempty(prepToPlot)
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        if isempty(prepIndices{r})
            prepIndices{r} = ...
                contains(PrepTypePerROI{r},prepToPlot);
            prepIndicesFly{r} = ...
                contains(PrepType{r},prepToPlot);
        else
            prepIndices{r} = (prepIndices{r} | ...
                contains(PrepTypePerROI{r},prepToPlot));
            prepIndicesFly{r} = (prepIndicesFly{r} | ...
                contains(PrepType{r},prepToPlot));
        end
    end
else
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        prepIndices{r} = true(size(PrepTypePerROI{r}));
        prepIndicesFly{r} = true(size(PrepType{r}));
    end
end
if ~isempty(gcampsToPlot)
    for p = 1:length(gcampsToPlot) %for each selected GCaMP type
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            if isempty(gcIndices{r})
                gcIndices{r} = ...
                    ismember(GCaMPTypePerROI{r},gcampsToPlot{p});
                gcIndicesFly{r} = ...
                    ismember(GCaMPType{r},gcampsToPlot{p});
            else
                gcIndices{r} = (gcIndices{r} | ...
                    ismember(GCaMPTypePerROI{r},gcampsToPlot{p}));
                gcIndicesFly{r} = (gcIndicesFly{r} | ...
                    ismember(GCaMPType{r},gcampsToPlot{p}));
            end
        end
    end
else
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        gcIndices{r} = true(size(GCaMPTypePerROI{r}));
        gcIndicesFly{r} = true(size(GCaMPType{r}));
    end
end

switch sideToPlot
    case 'Both'
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            sideIndices{r} = true(size(SomaSidePerROI{r}));
            sideIndicesFly{r} = true(size(SomaSide{r}));
        end
        switch normType
            case 'Zscore'
                dataToPlot = OFRes_1st_Z_R;
                aucToPlot = Auc_1st_bs_Z_R;
            case 'Max'
                dataToPlot = STResP_normM;
            case 'Raw'
                dataToPlot = OFRes_1st;
            case 'Znorm'
                dataToPlot = OFRes_1st_Znorm;
            otherwise
                error(['Normalization type not specified correctly. ',...
                    'Use normType=[''Max'',''Zscore'',''Znorm'', ',...
                    'or ''Raw'']'])
        end
    case {'L','R'}
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            if isempty(sideIndices{r})
                sideIndices{r} = ...
                    contains(SomaSidePerROI{r},sideToPlot);
                sideIndicesFly{r} = ...
                    contains(SomaSide{r},sideToPlot);
            else
                sideIndices{r} = (sideIndices{r} | ...
                    contains(SomaSidePerROI{r},sideToPlot));
                sideIndicesFly{r} = (sideIndicesFly{r} | ...
                    contains(SomaSide{r},sideToPlot));
            end
        end
        switch normType
            case 'Zscore'
                dataToPlot = OFRes_Raw_1st_Z;
            case 'Max'
                dataToPlot = OFRes_Raw_1st_normM;
            case 'Raw'
                dataToPlot1st = OFRes_Raw_1st;
            otherwise
                error(['Normalization type not specified correctly. ',...
                    'Use normType=[''Max'',''Zscore'', ',...
                    'or ''Raw'']'])
        end
    otherwise
        error(['Plotting side not specified correctly. ',...
            'Use sideToPlot = [''L'',''R'' or ''Both'']'])
end
% -------------------------------------------------------------------------
% Extended Data Figure 8a
f=plotSelectDFFCompound(dataToPlot,neuronsToPlot,...
    stimToPlot,minFrameRate,cReg,prepIndices,prepIndicesFly,...
    gcIndices,gcIndicesFly,sideIndices,sideIndicesFly,Stimuli,...
    regions,limsY,plotType);
% -------------------------------------------------------------------------
% Extended Data Figure 8b
stimToPlot=find(contains(Stimuli,{'Trans0-71','Trans1-42','Trans2-14'})&...
    contains(Stimuli,{'5px'}));
stimToPlot(contains(Stimuli(stimToPlot),{'Rot0','Rot1_CW'})) = [];
neuronsToPlot = [1,2];

f=plotSelectDFFCompound(dataToPlot,neuronsToPlot,...
    stimToPlot,minFrameRate,cReg,prepIndices,prepIndicesFly,...
    gcIndices,gcIndicesFly,sideIndices,sideIndicesFly,Stimuli,...
    regions,limsY,plotType);
% -------------------------------------------------------------------------

%% Extended Data Figure 8c

% which neurons to plot. Defined by their index in the regions cell array.
% Plotted in the same order. Leave blank to plot all neurons.

neuronsToPlot = [6,7];

% Which calcium indicator(s) to plot. Defined by the genotype field in OF
% struct array. Leave blank to plot everything. Note that it is possible to
% pool data from multiple indicators (e.g. gcampsToPlot = {GC6f','GC7f'})
% Accepted arguments: {'GC6f','GC7f','sytGC7f'} Default: gcampsToPlot = {};
gcampsToPlot = {'GC6f'};

% Which prep types to plot. Defined by the prepType field in OF struct
% array. Leave blank to plot all experiments.
% Accepted arguments: ['walking','immobilized'] Default: prepToPlot = '';
prepToPlot = '';

% Neurons from which hemisphere to plot. Defined by the somaSide field in
% OF struct array. Type 'Both' to plot all neurons. In that case, responses
% from Right hemisphere will be flipped to match the Left hemisphere neural
% responses and plotted together.
% Accepted arguments: ['L','R','Both'] Default: sideToPlot = 'Both';
sideToPlot = 'Both';

% specify normalization type. Type 'Raw' for no normalization.
% 'Max' will normalize the data per ROI to the maximum response recorded in
% that ROI (see variable maxes).
% 'Zscore' will zscore the data per ROI based on the mean and standard
% deviation of the entire recording session done on that ROI (see variables
% mu and sigma). 'Znorm' normalizes Z-score data by the max Z-score
% observed per fly per neuron
% Accepted arguments: ['Max','Zscore','Raw','Znorm'] Default: 'Raw';
normType = 'Zscore';
%----------------------------end of user input-----------------------------

clearvars l lh h stats

if isempty(stimToPlot)
    stimToPlot = 1:length(Stimuli);
end
if isempty(neuronsToPlot)
    neuronsToPlot = 1:length(regions);
end

prepIndices = cell(length(regions),1);
gcIndices = cell(length(regions),1);
sideIndices = cell(length(regions),1);
prepIndicesFly = cell(length(regions),1);
gcIndicesFly = cell(length(regions),1);
sideIndicesFly = cell(length(regions),1);

if ~isempty(prepToPlot)
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        if isempty(prepIndices{r})
            prepIndices{r} = ...
                contains(PrepTypePerROI{r},prepToPlot);
            prepIndicesFly{r} = ...
                contains(PrepType{r},prepToPlot);
        else
            prepIndices{r} = (prepIndices{r} | ...
                contains(PrepTypePerROI{r},prepToPlot));
            prepIndicesFly{r} = (prepIndicesFly{r} | ...
                contains(PrepType{r},prepToPlot));
        end
    end
else
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        prepIndices{r} = true(size(PrepTypePerROI{r}));
        prepIndicesFly{r} = true(size(PrepType{r}));
    end
end
if ~isempty(gcampsToPlot)
    for p = 1:length(gcampsToPlot) %for each selected GCaMP type
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            if isempty(gcIndices{r})
                gcIndices{r} = ...
                    ismember(GCaMPTypePerROI{r},gcampsToPlot{p});
                gcIndicesFly{r} = ...
                    ismember(GCaMPType{r},gcampsToPlot{p});
            else
                gcIndices{r} = (gcIndices{r} | ...
                    ismember(GCaMPTypePerROI{r},gcampsToPlot{p}));
                gcIndicesFly{r} = (gcIndicesFly{r} | ...
                    ismember(GCaMPType{r},gcampsToPlot{p}));
            end
        end
    end
else
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        gcIndices{r} = true(size(GCaMPTypePerROI{r}));
        gcIndicesFly{r} = true(size(GCaMPType{r}));
    end
end
switch sideToPlot
    case 'Both'
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            sideIndices{r} = true(size(SomaSidePerROI{r}));
            sideIndicesFly{r} = true(size(SomaSide{r}));
        end
        switch normType
            case 'Zscore'
                aucToPlot = Auc_1st_bs_Z_R;
                dataToPlot = OFRes_1st_Z_R;
            case 'Max'
                dataToPlot = Auc_1st_normM;
            case 'Raw'
                dataToPlot1st = Auc_1st_bs;
                dataToPlot2nd = Auc_2nd_bs;
            otherwise
                error(['Normalization type not specified correctly. ',...
                    'Use normType=[''Max'',''Zscore'', ',...
                    'or ''Raw'']'])
        end
    case {'L','R'}
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            if isempty(sideIndices{r})
                sideIndices{r} = ...
                    contains(SomaSidePerROI{r},sideToPlot);
                sideIndicesFly{r} = ...
                    contains(SomaSide{r},sideToPlot);
            else
                sideIndices{r} = (sideIndices{r} | ...
                    contains(SomaSidePerROI{r},sideToPlot));
                sideIndicesFly{r} = (sideIndicesFly{r} | ...
                    contains(SomaSide{r},sideToPlot));
            end
        end
        switch normType
            case 'Zscore'
                dataToPlot1st = Auc_1st_bs_Z;
                dataToPlot2nd = Auc_2nd_bs_Z;
            case 'Max'
                dataToPlot1st = Auc_1st_normM;
                dataToPlot2nd = Auc_2nd_normM;
            case 'Raw'
                dataToPlot1st = Auc_1st_bs;
                dataToPlot2nd = Auc_2nd_bs;
            otherwise
                error(['Normalization type not specified correctly. ',...
                    'Use normType=[''Max'',''Zscore'', ',...
                    'or ''Raw'']'])
        end
    otherwise
        error(['Plotting side not specified correctly. ',...
            'Use sideToPlot = [''L'',''R'' or ''Both'']'])
end

StimXLabel = {'zero','low','int.','fast'};
limsY = [0 1.2];

sp=0.5;
f=figure('WindowState','maximized');
x = zeros(size(neuronsToPlot,2),size(stimToPlotCW,2));
DSI = cell(size(neuronsToPlot,2),size(stimToPlotCW,2));
DSR = cell(size(neuronsToPlot,2),size(stimToPlotCW,2));
DSRnormM = cell(size(neuronsToPlot,2),size(stimToPlotCW,1));
DSRnormMA = cell(size(neuronsToPlot,2),size(stimToPlotCW,1));
my = zeros(1,size(neuronsToPlot,2));

for n = 1:length(neuronsToPlot)
    s = neuronsToPlot(n);
    iPlot = prepIndices{s} & gcIndices{s} & sideIndices{s};
    for k = 1:length(stimToPlotCW)
        mCW = stimToPlotCW(k);
        mCCW = stimToPlotCCW(k);
        x(n,k) = (n-1)*length(stimToPlotCW) + 1 + (k-1)*sp;
        

        StimXLabel{k} = Stimuli{stimToPlotCW(k)}(1:end-9);

        switch regions{s}
            case {'H2axon'}
                rPD = mean(dataToPlot{s,mCCW}(iPeak,iPlot));
                rNPD = mean(dataToPlot{s,mCW}(iPeak,iPlot));
            otherwise
                rPD = mean(dataToPlot{s,mCW}(iPeak,iPlot));
                rNPD = mean(dataToPlot{s,mCCW}(iPeak,iPlot));
        end
        DSI{n,k}=(rPD-rNPD)./(abs(rPD)+abs(rNPD));
        DSR{n,k}=(rPD-rNPD);
        % save the maximum y values for significance plot alignment
        maxbl = nanmax(DSR{n,k});
        if isempty(maxbl)
            maxbl = 0;
        end
        my(n) = max(my(n),maxbl);
    end
    maxnorm = zeros(1,size(DSR{n,k},2));
    for i = 1:size(DSR{n,k},2)
        for k = 1:length(stimToPlotCW)
            if i <= size(DSR{n,k},2)
                maxnorm(1,i) = max([maxnorm(1,i), abs(DSR{n,k}(i))]);
            end
        end
    end
    for k = 1:length(stimToPlotCW)
        DSRnormM{n,k}=DSR{n,k}./maxnorm(1,1:size(DSR{n,k},2));
        DSRnormMA{n,k}=abs(DSR{n,k})./maxnorm(1,1:size(DSR{n,k},2));
        DSI{n,k}(isnan(DSI{n,k})) = [];
        DSR{n,k}(isnan(DSR{n,k})) = [];
        DSRnormM{n,k}(isnan(DSRnormM{n,k})) = [];
        DSRnormMA{n,k}(isnan(DSRnormMA{n,k})) = [];
    end
    for k = 1:4
        sti = k+5;
        x(n,k) = (n-1)*length(stimToPlot) + 1 + (k-1)*sp;
        h(n,k) = notBoxPlot(DSRnormMA{n,sti},x(n,k));hold on
        set(h(n,k).data,'MarkerFaceColor',cReg(s,:))
    end
end
limsX = xlim;
set(gca,'XTick',unique(x),'XTickLabel',StimXLabel,...
    'FontSize',14,'TickLabelInterpreter','none',...
    'XTickLabelRotation',45)
set(gca,'YLim',limsY)
ylabel('Asymmetry Index')

for m=1:length(neuronsToPlot)
    hold on
    l(m)=plot(nan,nan,'Color',cReg(neuronsToPlot(m),:),'LineWidth',4);
end
set(gca,'XLim',limsX)
lh=legend(l,regions(neuronsToPlot),'Box','off','FontSize',20,...
    'Location','Best','Interpreter','none');

%% Extended Data Figure 8d

% which neurons to plot. Defined by their index in the regions cell array.
% Plotted in the same order. Leave blank to plot all neurons.

neuronsToPlot = [1,2];

% Which calcium indicator(s) to plot. Defined by the genotype field in OF
% struct array. Leave blank to plot everything. Note that it is possible to
% pool data from multiple indicators (e.g. gcampsToPlot = {GC6f','GC7f'})
% Accepted arguments: {'GC6f','GC7f','sytGC7f'} Default: gcampsToPlot = {};
gcampsToPlot = {'GC6f'};

% Which prep types to plot. Defined by the prepType field in OF struct
% array. Leave blank to plot all experiments.
% Accepted arguments: ['walking','immobilized'] Default: prepToPlot = '';
prepToPlot = '';

% Neurons from which hemisphere to plot. Defined by the somaSide field in
% OF struct array. Type 'Both' to plot all neurons. In that case, responses
% from Right hemisphere will be flipped to match the Left hemisphere neural
% responses and plotted together.
% Accepted arguments: ['L','R','Both'] Default: sideToPlot = 'Both';
sideToPlot = 'Both';

% specify normalization type. Type 'Raw' for no normalization.
% 'Max' will normalize the data per ROI to the maximum response recorded in
% that ROI (see variable maxes).
% 'Zscore' will zscore the data per ROI based on the mean and standard
% deviation of the entire recording session done on that ROI (see variables
% mu and sigma). 'Znorm' normalizes Z-score data by the max Z-score
% observed per fly per neuron
% Accepted arguments: ['Max','Zscore','Raw','Znorm'] Default: 'Raw';
normType = 'Zscore';
%----------------------------end of user input-----------------------------

clearvars l lh h stats

if isempty(stimToPlot)
    stimToPlot = 1:length(Stimuli);
end
if isempty(neuronsToPlot)
    neuronsToPlot = 1:length(regions);
end

prepIndices = cell(length(regions),1);
gcIndices = cell(length(regions),1);
sideIndices = cell(length(regions),1);
prepIndicesFly = cell(length(regions),1);
gcIndicesFly = cell(length(regions),1);
sideIndicesFly = cell(length(regions),1);

if ~isempty(prepToPlot)
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        if isempty(prepIndices{r})
            prepIndices{r} = ...
                contains(PrepTypePerROI{r},prepToPlot);
            prepIndicesFly{r} = ...
                contains(PrepType{r},prepToPlot);
        else
            prepIndices{r} = (prepIndices{r} | ...
                contains(PrepTypePerROI{r},prepToPlot));
            prepIndicesFly{r} = (prepIndicesFly{r} | ...
                contains(PrepType{r},prepToPlot));
        end
    end
else
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        prepIndices{r} = true(size(PrepTypePerROI{r}));
        prepIndicesFly{r} = true(size(PrepType{r}));
    end
end
if ~isempty(gcampsToPlot)
    for p = 1:length(gcampsToPlot) %for each selected GCaMP type
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            if isempty(gcIndices{r})
                gcIndices{r} = ...
                    ismember(GCaMPTypePerROI{r},gcampsToPlot{p});
                gcIndicesFly{r} = ...
                    ismember(GCaMPType{r},gcampsToPlot{p});
            else
                gcIndices{r} = (gcIndices{r} | ...
                    ismember(GCaMPTypePerROI{r},gcampsToPlot{p}));
                gcIndicesFly{r} = (gcIndicesFly{r} | ...
                    ismember(GCaMPType{r},gcampsToPlot{p}));
            end
        end
    end
else
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        gcIndices{r} = true(size(GCaMPTypePerROI{r}));
        gcIndicesFly{r} = true(size(GCaMPType{r}));
    end
end
switch sideToPlot
    case 'Both'
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            sideIndices{r} = true(size(SomaSidePerROI{r}));
            sideIndicesFly{r} = true(size(SomaSide{r}));
        end
        switch normType
            case 'Zscore'
                aucToPlot = Auc_1st_bs_Z_R;
                dataToPlot = OFRes_1st_Z_R;
            case 'Max'
                dataToPlot = Auc_1st_normM;
            case 'Raw'
                dataToPlot1st = Auc_1st_bs;
                dataToPlot2nd = Auc_2nd_bs;
            otherwise
                error(['Normalization type not specified correctly. ',...
                    'Use normType=[''Max'',''Zscore'', ',...
                    'or ''Raw'']'])
        end
    case {'L','R'}
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            if isempty(sideIndices{r})
                sideIndices{r} = ...
                    contains(SomaSidePerROI{r},sideToPlot);
                sideIndicesFly{r} = ...
                    contains(SomaSide{r},sideToPlot);
            else
                sideIndices{r} = (sideIndices{r} | ...
                    contains(SomaSidePerROI{r},sideToPlot));
                sideIndicesFly{r} = (sideIndicesFly{r} | ...
                    contains(SomaSide{r},sideToPlot));
            end
        end
        switch normType
            case 'Zscore'
                dataToPlot1st = Auc_1st_bs_Z;
                dataToPlot2nd = Auc_2nd_bs_Z;
            case 'Max'
                dataToPlot1st = Auc_1st_normM;
                dataToPlot2nd = Auc_2nd_normM;
            case 'Raw'
                dataToPlot1st = Auc_1st_bs;
                dataToPlot2nd = Auc_2nd_bs;
            otherwise
                error(['Normalization type not specified correctly. ',...
                    'Use normType=[''Max'',''Zscore'', ',...
                    'or ''Raw'']'])
        end
    otherwise
        error(['Plotting side not specified correctly. ',...
            'Use sideToPlot = [''L'',''R'' or ''Both'']'])
end

StimXLabel = {'zero','low','int.','fast'};
limsY = [0 1.2];

sp=0.5;
f=figure('WindowState','maximized');
x = zeros(size(neuronsToPlot,2),size(stimToPlotCW,2));
DSI = cell(size(neuronsToPlot,2),size(stimToPlotCW,2));
DSR = cell(size(neuronsToPlot,2),size(stimToPlotCW,2));
DSRnormM = cell(size(neuronsToPlot,2),size(stimToPlotCW,1));
DSRnormMA = cell(size(neuronsToPlot,2),size(stimToPlotCW,1));
my = zeros(1,size(neuronsToPlot,2));

for n = 1:length(neuronsToPlot)
    s = neuronsToPlot(n);
    iPlot = prepIndices{s} & gcIndices{s} & sideIndices{s};
    for k = 1:length(stimToPlotCW)
        mCW = stimToPlotCW(k);
        mCCW = stimToPlotCCW(k);
        x(n,k) = (n-1)*length(stimToPlotCW) + 1 + (k-1)*sp;
        

        StimXLabel{k} = Stimuli{stimToPlotCW(k)}(1:end-9);

        switch regions{s}
            case {'H2axon'}
                rPD = mean(dataToPlot{s,mCCW}(iPeak,iPlot));
                rNPD = mean(dataToPlot{s,mCW}(iPeak,iPlot));
            otherwise
                rPD = mean(dataToPlot{s,mCW}(iPeak,iPlot));
                rNPD = mean(dataToPlot{s,mCCW}(iPeak,iPlot));
        end
        DSI{n,k}=(rPD-rNPD)./(abs(rPD)+abs(rNPD));
        DSR{n,k}=(rPD-rNPD);
        % save the maximum y values for significance plot alignment
        maxbl = nanmax(DSR{n,k});
        if isempty(maxbl)
            maxbl = 0;
        end
        my(n) = max(my(n),maxbl);
    end
    maxnorm = zeros(1,size(DSR{n,k},2));
    for i = 1:size(DSR{n,k},2)
        for k = 1:length(stimToPlotCW)
            if i <= size(DSR{n,k},2)
                maxnorm(1,i) = max([maxnorm(1,i), abs(DSR{n,k}(i))]);
            end
        end
    end
    for k = 1:length(stimToPlotCW)
        DSRnormM{n,k}=DSR{n,k}./maxnorm(1,1:size(DSR{n,k},2));
        DSRnormMA{n,k}=abs(DSR{n,k})./maxnorm(1,1:size(DSR{n,k},2));
        DSI{n,k}(isnan(DSI{n,k})) = [];
        DSR{n,k}(isnan(DSR{n,k})) = [];
        DSRnormM{n,k}(isnan(DSRnormM{n,k})) = [];
        DSRnormMA{n,k}(isnan(DSRnormMA{n,k})) = [];
    end
    for k = 1:4
        sti = k+5;
        x(n,k) = (n-1)*length(stimToPlot) + 1 + (k-1)*sp;
        h(n,k) = notBoxPlot(DSRnormMA{n,sti},x(n,k));hold on
        set(h(n,k).data,'MarkerFaceColor',cReg(s,:))
    end
end
limsX = xlim;
set(gca,'XTick',unique(x),'XTickLabel',StimXLabel,...
    'FontSize',14,'TickLabelInterpreter','none',...
    'XTickLabelRotation',45)
set(gca,'YLim',limsY)
ylabel('Asymmetry Index')

for m=1:length(neuronsToPlot)
    hold on
    l(m)=plot(nan,nan,'Color',cReg(neuronsToPlot(m),:),'LineWidth',4);
end
set(gca,'XLim',limsX)
lh=legend(l,regions(neuronsToPlot),'Box','off','FontSize',20,...
    'Location','Best','Interpreter','none');

%% Extended Data Figure 8e


% Which compound OF stimuli to plot. Defined by their index in the Stimuli
% cell array. Plotted in the same order. Leave blank to plot all stimuli.

% stimToPlot = find(contains(Stimuli,{'Trans1-42Rot1'})&contains(Stimuli,{'5px'}));
stimToPlot = find(contains(Stimuli,{'Flicker'})&contains(Stimuli,{'5px'}));

% which neurons to plot. Defined by their index in the regions cell array.
% Plotted in the same order. Leave blank to plot all neurons.

neuronsToPlot = [6,7,9,10,1,2];

% Which calcium indicator(s) to plot. Defined by the genotype field in OF
% struct array. Leave blank to plot everything. Note that it is possible to
% pool data from multiple indicators (e.g. gcampsToPlot = {GC6f','GC7f'})
% Accepted arguments: {'GC6f','GC7f','sytGC7f'} Default: gcampsToPlot = {};
gcampsToPlot = {'GC6f'};

% Which prep types to plot. Defined by the prepType field in OF struct
% array. Leave blank to plot all experiments.
% Accepted arguments: ['walking','immobilized'] Default: prepToPlot = '';
prepToPlot = '';

% Neurons from which hemisphere to plot. Defined by the somaSide field in
% OF struct array. Type 'Both' to plot all neurons. In that case, responses
% from Right hemisphere will be flipped to match the Left hemisphere neural
% responses and plotted together.
% Accepted arguments: ['L','R','Both'] Default: sideToPlot = 'Both';
sideToPlot = 'Both';

% specify normalization type. Type 'Raw' for no normalization.
% 'Max' will normalize the data per ROI to the maximum response recorded in
% that ROI (see variable maxes).
% 'Zscore' will zscore the data per ROI based on the mean and standard
% deviation of the entire recording session done on that ROI (see variables
% mu and sigma). 'Znorm' normalizes Z-score data by the max Z-score
% observed per fly per neuron
% Accepted arguments: ['Max','Zscore','Raw','Znorm'] Default: 'Raw';
normType = 'Zscore';

%----------------------------end of user input-----------------------------

peakbin = 500; % time bin for the averaged peak response (ms)
clearvars l lh h stats

% Define the subset of data that will be plotted. This will be used in the
% plotSelectAUC function to plot only the specified subset of data as well
% as to calculate and display number of flies & ROIs being plotted
if isempty(stimToPlot)
    stimToPlot = 1:length(Stimuli);
end
if isempty(neuronsToPlot)
    neuronsToPlot = 1:length(regions);
end

prepIndices = cell(length(regions),1);
gcIndices = cell(length(regions),1);
sideIndices = cell(length(regions),1);
prepIndicesFly = cell(length(regions),1);
gcIndicesFly = cell(length(regions),1);
sideIndicesFly = cell(length(regions),1);

if ~isempty(prepToPlot)
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        if isempty(prepIndices{r})
            prepIndices{r} = ...
                contains(PrepTypePerROI{r},prepToPlot);
            prepIndicesFly{r} = ...
                contains(PrepType{r},prepToPlot);
        else
            prepIndices{r} = (prepIndices{r} | ...
                contains(PrepTypePerROI{r},prepToPlot));
            prepIndicesFly{r} = (prepIndicesFly{r} | ...
                contains(PrepType{r},prepToPlot));
        end
    end
else
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        prepIndices{r} = true(size(PrepTypePerROI{r}));
        prepIndicesFly{r} = true(size(PrepType{r}));
    end
end
if ~isempty(gcampsToPlot)
    for p = 1:length(gcampsToPlot) %for each selected GCaMP type
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            if isempty(gcIndices{r})
                gcIndices{r} = ...
                    ismember(GCaMPTypePerROI{r},gcampsToPlot{p});
                gcIndicesFly{r} = ...
                    ismember(GCaMPType{r},gcampsToPlot{p});
            else
                gcIndices{r} = (gcIndices{r} | ...
                    ismember(GCaMPTypePerROI{r},gcampsToPlot{p}));
                gcIndicesFly{r} = (gcIndicesFly{r} | ...
                    ismember(GCaMPType{r},gcampsToPlot{p}));
            end
        end
    end
else
    for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
        gcIndices{r} = true(size(GCaMPTypePerROI{r}));
        gcIndicesFly{r} = true(size(GCaMPType{r}));
    end
end
switch sideToPlot
    case 'Both'
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            sideIndices{r} = true(size(SomaSidePerROI{r}));
            sideIndicesFly{r} = true(size(SomaSide{r}));
        end
        switch normType
            case 'Zscore'
                aucToPlot = Auc_1st_bs_Z_R;
                dataToPlot = OFRes_1st_Z_R;
            case 'Max'
                dataToPlot = Auc_1st_normM;
            case 'Raw'
                dataToPlot1st = Auc_1st_bs;
                dataToPlot2nd = Auc_2nd_bs;
            otherwise
                error(['Normalization type not specified correctly. ',...
                    'Use normType=[''Max'',''Zscore'', ',...
                    'or ''Raw'']'])
        end
    case {'L','R'}
        for r = neuronsToPlot(1:length(neuronsToPlot)) %for each region
            if isempty(sideIndices{r})
                sideIndices{r} = ...
                    contains(SomaSidePerROI{r},sideToPlot);
                sideIndicesFly{r} = ...
                    contains(SomaSide{r},sideToPlot);
            else
                sideIndices{r} = (sideIndices{r} | ...
                    contains(SomaSidePerROI{r},sideToPlot));
                sideIndicesFly{r} = (sideIndicesFly{r} | ...
                    contains(SomaSide{r},sideToPlot));
            end
        end
        switch normType
            case 'Zscore'
                dataToPlot1st = Auc_1st_bs_Z;
                dataToPlot2nd = Auc_2nd_bs_Z;
            case 'Max'
                dataToPlot1st = Auc_1st_normM;
                dataToPlot2nd = Auc_2nd_normM;
            case 'Raw'
                dataToPlot1st = Auc_1st_bs;
                dataToPlot2nd = Auc_2nd_bs;
            otherwise
                error(['Normalization type not specified correctly. ',...
                    'Use normType=[''Max'',''Zscore'', ',...
                    'or ''Raw'']'])
        end
    otherwise
        error(['Plotting side not specified correctly. ',...
            'Use sideToPlot = [''L'',''R'' or ''Both'']'])
end

frameNum = size(dataToPlot{1,1},1);
TimeAxisPD = linspace(0,frameNum/(minFrameRate/1000),frameNum);
limsY = [-5 30];
iPeak = find(TimeAxisPD>4-peakbin/1000,1,'first'):...
    find(TimeAxisPD<4,1,'last');

sp=0.5;
f=figure('WindowState','maximized');
StimXLabel = {'Left','Right'};
x = zeros(size(neuronsToPlot,2),size(stimToPlotCW,2));
DSI = cell(size(neuronsToPlot,2),size(stimToPlotCW,2));
DSR = cell(size(neuronsToPlot,2),size(stimToPlotCW,2));
DSRnormM = cell(size(neuronsToPlot,2),size(stimToPlotCW,1));
DSRnormMA = cell(size(neuronsToPlot,2),size(stimToPlotCW,1));
NPD = cell(size(neuronsToPlot,2),size(stimToPlotCW,2));
PD = cell(size(neuronsToPlot,2),size(stimToPlotCW,2));
my = zeros(1,size(neuronsToPlot,2));

for n = 1:length(neuronsToPlot)
    s = neuronsToPlot(n);
    iPlot = prepIndices{s} & gcIndices{s} & sideIndices{s};
    for k = 1:length(stimToPlotCW)
        mCW = stimToPlotCW(k);
        mCCW = stimToPlotCCW(k);
        switch regions{s}
            case {'H2axon'}
                rPD = mean(dataToPlot{s,mCCW}(iPeak,iPlot));
                rNPD = mean(dataToPlot{s,mCW}(iPeak,iPlot));
            otherwise
                rPD = mean(dataToPlot{s,mCW}(iPeak,iPlot));
                rNPD = mean(dataToPlot{s,mCCW}(iPeak,iPlot));
        end
        DSI{n,k}=(rPD-rNPD)./(abs(rPD)+abs(rNPD));
        DSR{n,k}=(rPD-rNPD);
        PD{n,k}=(rPD);
        NPD{n,k}=(rNPD);
        maxbl = nanmax(DSR{n,k});
        if isempty(maxbl)
            maxbl = 0;
        end
        my(n) = max(my(n),maxbl);
        NPD{n,k}(isnan(NPD{n,k})) = [];
        PD{n,k}(isnan(PD{n,k})) = [];
    end
    for k = 1:length(stimToPlot)
        if sum(stimToPlotCW==stimToPlot(k))
            sti = find(stimToPlotCW==stimToPlot(k));
            x(n,k) = (n-1)*length(stimToPlot) + 1 + (k-1)*sp;
            h(n,k) = notBoxPlot(NPD{n,sti},x(n,k));hold on
            set(h(n,k).data,'MarkerFaceColor',cReg(s,:))
        else
            sti = find(stimToPlotCCW==stimToPlot(k));
            x(n,k) = (n-1)*length(stimToPlot) + 1 + (k-1)*sp;
            h(n,k) = notBoxPlot(PD{n,sti},x(n,k));hold on
            set(h(n,k).data,'MarkerFaceColor',cReg(s,:))
        end
    end
end
limsX = xlim;
line(limsX,[0,0],'LineStyle','--','Color','k','LineWidth',1)
set(gca,'XTick',unique(x),'XTickLabel',StimXLabel,...
    'FontSize',14,'TickLabelInterpreter','none',...
    'XTickLabelRotation',45)
set(gca,'YLim',limsY)
ylabel('Deviated translation Response (z-score)')
title('flicker')
for m=1:length(neuronsToPlot)
    hold on
    l(m)=plot(nan,nan,'Color',cReg(neuronsToPlot(m),:),'LineWidth',4);
end
set(gca,'XLim',limsX)
lh=legend(l,regions(neuronsToPlot),'Box','off','FontSize',20,...
    'Location','Best','Interpreter','none');

% -------------------------------------------------------------------------
% zero speed
f=figure('WindowState','maximized');
stimToPlot = find(contains(Stimuli,{'Trans0Rot1'})&contains(Stimuli,{'5px'}));
x = zeros(size(neuronsToPlot,2),size(stimToPlotCW,2));
DSI = cell(size(neuronsToPlot,2),size(stimToPlotCW,2));
DSR = cell(size(neuronsToPlot,2),size(stimToPlotCW,2));
DSRnormM = cell(size(neuronsToPlot,2),size(stimToPlotCW,1));
DSRnormMA = cell(size(neuronsToPlot,2),size(stimToPlotCW,1));
NPD = cell(size(neuronsToPlot,2),size(stimToPlotCW,2));
PD = cell(size(neuronsToPlot,2),size(stimToPlotCW,2));
my = zeros(1,size(neuronsToPlot,2));

for n = 1:length(neuronsToPlot)
    s = neuronsToPlot(n);
    iPlot = prepIndices{s} & gcIndices{s} & sideIndices{s};
    for k = 1:length(stimToPlotCW)
        mCW = stimToPlotCW(k);
        mCCW = stimToPlotCCW(k);
        switch regions{s}
            case {'H2axon'}
                rPD = mean(dataToPlot{s,mCCW}(iPeak,iPlot));
                rNPD = mean(dataToPlot{s,mCW}(iPeak,iPlot));
            otherwise
                rPD = mean(dataToPlot{s,mCW}(iPeak,iPlot));
                rNPD = mean(dataToPlot{s,mCCW}(iPeak,iPlot));
        end
        DSI{n,k}=(rPD-rNPD)./(abs(rPD)+abs(rNPD));
        DSR{n,k}=(rPD-rNPD);
        PD{n,k}=(rPD);
        NPD{n,k}=(rNPD);
        maxbl = nanmax(DSR{n,k});
        if isempty(maxbl)
            maxbl = 0;
        end
        my(n) = max(my(n),maxbl);
        NPD{n,k}(isnan(NPD{n,k})) = [];
        PD{n,k}(isnan(PD{n,k})) = [];
    end
    for k = 1:length(stimToPlot)
        if sum(stimToPlotCW==stimToPlot(k))
            sti = find(stimToPlotCW==stimToPlot(k));
            x(n,k) = (n-1)*length(stimToPlot) + 1 + (k-1)*sp;
            h(n,k) = notBoxPlot(NPD{n,sti},x(n,k));hold on
            set(h(n,k).data,'MarkerFaceColor',cReg(s,:))
        else
            sti = find(stimToPlotCCW==stimToPlot(k));
            x(n,k) = (n-1)*length(stimToPlot) + 1 + (k-1)*sp;
            h(n,k) = notBoxPlot(PD{n,sti},x(n,k));hold on
            set(h(n,k).data,'MarkerFaceColor',cReg(s,:))
        end
    end
end
limsX = xlim;
line(limsX,[0,0],'LineStyle','--','Color','k','LineWidth',1)
set(gca,'XTick',unique(x),'XTickLabel',StimXLabel,...
    'FontSize',14,'TickLabelInterpreter','none',...
    'XTickLabelRotation',45)
set(gca,'YLim',limsY)
ylabel('Deviated translation Response (z-score)')
title('zero speed')
for m=1:length(neuronsToPlot)
    hold on
    l(m)=plot(nan,nan,'Color',cReg(neuronsToPlot(m),:),'LineWidth',4);
end
set(gca,'XLim',limsX)
lh=legend(l,regions(neuronsToPlot),'Box','off','FontSize',20,...
    'Location','Best','Interpreter','none');
% -------------------------------------------------------------------------
% low speed
f=figure('WindowState','maximized');
stimToPlot = find(contains(Stimuli,{'Trans0-71Rot1'})&contains(Stimuli,{'5px'}));
x = zeros(size(neuronsToPlot,2),size(stimToPlotCW,2));
DSI = cell(size(neuronsToPlot,2),size(stimToPlotCW,2));
DSR = cell(size(neuronsToPlot,2),size(stimToPlotCW,2));
DSRnormM = cell(size(neuronsToPlot,2),size(stimToPlotCW,1));
DSRnormMA = cell(size(neuronsToPlot,2),size(stimToPlotCW,1));
NPD = cell(size(neuronsToPlot,2),size(stimToPlotCW,2));
PD = cell(size(neuronsToPlot,2),size(stimToPlotCW,2));
my = zeros(1,size(neuronsToPlot,2));

for n = 1:length(neuronsToPlot)
    s = neuronsToPlot(n);
    iPlot = prepIndices{s} & gcIndices{s} & sideIndices{s};
    for k = 1:length(stimToPlotCW)
        mCW = stimToPlotCW(k);
        mCCW = stimToPlotCCW(k);
        switch regions{s}
            case {'H2axon'}
                rPD = mean(dataToPlot{s,mCCW}(iPeak,iPlot));
                rNPD = mean(dataToPlot{s,mCW}(iPeak,iPlot));
            otherwise
                rPD = mean(dataToPlot{s,mCW}(iPeak,iPlot));
                rNPD = mean(dataToPlot{s,mCCW}(iPeak,iPlot));
        end
        DSI{n,k}=(rPD-rNPD)./(abs(rPD)+abs(rNPD));
        DSR{n,k}=(rPD-rNPD);
        PD{n,k}=(rPD);
        NPD{n,k}=(rNPD);
        maxbl = nanmax(DSR{n,k});
        if isempty(maxbl)
            maxbl = 0;
        end
        my(n) = max(my(n),maxbl);
        NPD{n,k}(isnan(NPD{n,k})) = [];
        PD{n,k}(isnan(PD{n,k})) = [];
    end
    for k = 1:length(stimToPlot)
        if sum(stimToPlotCW==stimToPlot(k))
            sti = find(stimToPlotCW==stimToPlot(k));
            x(n,k) = (n-1)*length(stimToPlot) + 1 + (k-1)*sp;
            h(n,k) = notBoxPlot(NPD{n,sti},x(n,k));hold on
            set(h(n,k).data,'MarkerFaceColor',cReg(s,:))
        else
            sti = find(stimToPlotCCW==stimToPlot(k));
            x(n,k) = (n-1)*length(stimToPlot) + 1 + (k-1)*sp;
            h(n,k) = notBoxPlot(PD{n,sti},x(n,k));hold on
            set(h(n,k).data,'MarkerFaceColor',cReg(s,:))
        end
    end
end
limsX = xlim;
line(limsX,[0,0],'LineStyle','--','Color','k','LineWidth',1)
set(gca,'XTick',unique(x),'XTickLabel',StimXLabel,...
    'FontSize',14,'TickLabelInterpreter','none',...
    'XTickLabelRotation',45)
set(gca,'YLim',limsY)
ylabel('Deviated translation Response (z-score)')
title('low speed')
for m=1:length(neuronsToPlot)
    hold on
    l(m)=plot(nan,nan,'Color',cReg(neuronsToPlot(m),:),'LineWidth',4);
end
set(gca,'XLim',limsX)
lh=legend(l,regions(neuronsToPlot),'Box','off','FontSize',20,...
    'Location','Best','Interpreter','none');
% -------------------------------------------------------------------------
% fast speed
f=figure('WindowState','maximized');
stimToPlot = find(contains(Stimuli,{'Trans2-14Rot1'})&contains(Stimuli,{'5px'}));
x = zeros(size(neuronsToPlot,2),size(stimToPlotCW,2));
DSI = cell(size(neuronsToPlot,2),size(stimToPlotCW,2));
DSR = cell(size(neuronsToPlot,2),size(stimToPlotCW,2));
DSRnormM = cell(size(neuronsToPlot,2),size(stimToPlotCW,1));
DSRnormMA = cell(size(neuronsToPlot,2),size(stimToPlotCW,1));
NPD = cell(size(neuronsToPlot,2),size(stimToPlotCW,2));
PD = cell(size(neuronsToPlot,2),size(stimToPlotCW,2));
my = zeros(1,size(neuronsToPlot,2));

for n = 1:length(neuronsToPlot)
    s = neuronsToPlot(n);
    iPlot = prepIndices{s} & gcIndices{s} & sideIndices{s};
    for k = 1:length(stimToPlotCW)
        mCW = stimToPlotCW(k);
        mCCW = stimToPlotCCW(k);
        switch regions{s}
            case {'H2axon'}
                rPD = mean(dataToPlot{s,mCCW}(iPeak,iPlot));
                rNPD = mean(dataToPlot{s,mCW}(iPeak,iPlot));
            otherwise
                rPD = mean(dataToPlot{s,mCW}(iPeak,iPlot));
                rNPD = mean(dataToPlot{s,mCCW}(iPeak,iPlot));
        end
        DSI{n,k}=(rPD-rNPD)./(abs(rPD)+abs(rNPD));
        DSR{n,k}=(rPD-rNPD);
        PD{n,k}=(rPD);
        NPD{n,k}=(rNPD);
        maxbl = nanmax(DSR{n,k});
        if isempty(maxbl)
            maxbl = 0;
        end
        my(n) = max(my(n),maxbl);
        NPD{n,k}(isnan(NPD{n,k})) = [];
        PD{n,k}(isnan(PD{n,k})) = [];
    end
    for k = 1:length(stimToPlot)
        if sum(stimToPlotCW==stimToPlot(k))
            sti = find(stimToPlotCW==stimToPlot(k));
            x(n,k) = (n-1)*length(stimToPlot) + 1 + (k-1)*sp;
            h(n,k) = notBoxPlot(NPD{n,sti},x(n,k));hold on
            set(h(n,k).data,'MarkerFaceColor',cReg(s,:))
        else
            sti = find(stimToPlotCCW==stimToPlot(k));
            x(n,k) = (n-1)*length(stimToPlot) + 1 + (k-1)*sp;
            h(n,k) = notBoxPlot(PD{n,sti},x(n,k));hold on
            set(h(n,k).data,'MarkerFaceColor',cReg(s,:))
        end
    end
end
limsX = xlim;
line(limsX,[0,0],'LineStyle','--','Color','k','LineWidth',1)
set(gca,'XTick',unique(x),'XTickLabel',StimXLabel,...
    'FontSize',14,'TickLabelInterpreter','none',...
    'XTickLabelRotation',45)
set(gca,'YLim',limsY)
ylabel('Deviated translation Response (z-score)')
title('fast speed')
for m=1:length(neuronsToPlot)
    hold on
    l(m)=plot(nan,nan,'Color',cReg(neuronsToPlot(m),:),'LineWidth',4);
end
set(gca,'XLim',limsX)
lh=legend(l,regions(neuronsToPlot),'Box','off','FontSize',20,...
    'Location','Best','Interpreter','none');
% -------------------------------------------------------------------------

% Optional Savefigures  
if SaveFigures
    % Create a timestamp string
    ts = datestr(now, 'yymmdd_HHMM');

    % Get all open figure handles
    figHandles = findall(0, 'Type', 'figure');

    % Sort them by ascending figure number so figure #1 is handled first, etc.
    [~, sortIdx] = sort([figHandles.Number]);
    figHandles = figHandles(sortIdx);

    % Loop through each figure handle in ascending numeric order
    for i = 1:length(figHandles)
        figure(figHandles(i));    % Make sure this figure is active
        drawnow;                  % Force an update of the figure

        % Fill the entire screen
        set(figHandles(i), 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

        drawnow; pause(1);        % Give time to resize

        % Build a filename that includes the figure number and timestamp
        figNum  = figHandles(i).Number;
        figName = sprintf('Figure_%d_%s', figNum, ts);

        % Save in PNG format
        saveas(figHandles(i), fullfile(saveDir, [figName, '.png']));

        % Optionally close
        close(figHandles(i));
    end
end
