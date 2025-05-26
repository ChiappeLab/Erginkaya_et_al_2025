function [SPKVR, SPKVRRand, SPKVRDist,SPKVF, pcav, pTypes, dtt] = GetSpikePCAMert(path, params)

flies = dir(path);
flies = flies(3:end);
dt = load([path flies(1).name '\DataLowRes.mat'], 'Flies');
seq = dt.Flies.Seq;
pTypes = unique(seq);
disp(path)
SPKVR = cell(length(flies),length(pTypes));
SPKVF = cell(length(flies),length(pTypes));
SPKVS = cell(length(flies),length(pTypes));
SPKVRRand = cell(length(flies),length(pTypes));
SPKVFRand = cell(length(flies),length(pTypes));
SPKVSRand = cell(length(flies),length(pTypes));
SPKVRDist = cell(length(flies),length(pTypes));
SPKTIMES = cell(length(flies),length(pTypes));
TOTALTIME = zeros(length(flies),length(seq));
VF = cell(length(flies),length(pTypes));
VR = cell(length(flies),length(pTypes));
VS = cell(length(flies),length(pTypes));
ACST = cell(length(flies),length(pTypes));
wi = 0; % wi = 20
dm = 1; % dm = 7
for n = 1 : length(flies)
    if(exist([path flies(n).name '\DataLowRes.mat'], 'file') == 2)
        dt = load([path flies(n).name '\DataLowRes.mat'], 'Flies');
        disp(['Running Fly ' num2str(n) ' of ' num2str(length(flies))])
        seq = dt.Flies.Seq;
        dt = dt.Flies.Data;
        totaltime = zeros(1,length(seq));
        for k = 1 : length(dt)
            for j = 1 : length(pTypes)
                spkvr = [];
                spkvf = [];
                spkvs = [];
                spkvrrand = [];
                spkvfrand = [];
                spkvsrand = [];
                spkdist = [];
                spikeTimes = [];
                Vf = []; Vs = []; Vr = []; acSt = [];
                switch seq{k}
                    case pTypes{j}
                        Vf = dt{k}.Vf;
                        Vs = dt{k}.Vs;
                        Vr = dt{k}.Vr;
                        wd = dt{k}.WallDist;
                        acSt = dt{k}.actState;
                        actst = acSt;
                        actst(wd < params.mDistWall) = 0;
                        [locs, ~, cmhSong, thr] = SortSpikes(Vr, actst, params);
                        [locsF, pksIF,pksEF, ~] = CullSpikesBasedOnVf(Vr, locs, cmhSong, thr, Vf, params);
                        [locsF, pksIF,pksEF, ~] = TemplateCompSpike(Vr,locsF,pksIF, pksEF, cmhSong, thr, params);
                        locsRand = randi([wi+dm+1 length(Vr)-(wi+dm+1)], length(locsF),1);
                        spikeTimes = [locsF-pksIF locsF locsF+pksEF];
                        spikeTimes(:,4) = nan;
                        totaltime(k) = length(Vf);
                        history = 30; % amount of frames preceding a saccade to determine if the fly was stationary or not
                        for t = 1:size(spikeTimes,1)
                            if spikeTimes(t,1)<= history
                                spikeTimes(t,4) = nan;
                            else
                                spikeTimes(t,4) = sum(acSt((spikeTimes(t,1)-history):(spikeTimes(t,1)-1)))/history;
                            end
%                             if spikeTimes(t,4) == 0 % Check the profile of saccades
%                                 if spikeTimes(t,3)+history >= length(Vr)
%                                     figure(1);
%                                     subplot(3,1,1);hold on;
%                                     plot(Vr((spikeTimes(t,1)-history):end))
%                                     subplot(3,1,2);hold on;
%                                     plot(Vf((spikeTimes(t,1)-history):end))
%                                     subplot(3,1,3);hold on;
%                                     plot(Vs((spikeTimes(t,1)-history):end))
%                                 else
%                                     figure(1);
%                                     subplot(3,1,1);hold on;
%                                     plot(Vr((spikeTimes(t,1)-history):(spikeTimes(t,3)+history)))
%                                     subplot(3,1,2);hold on;
%                                     plot(Vf((spikeTimes(t,1)-history):(spikeTimes(t,3)+history)))
%                                     subplot(3,1,3);hold on;
%                                     plot(Vs((spikeTimes(t,1)-history):(spikeTimes(t,3)+history)))
%                                 end
%                             end
                        end
                        if(length(locsF) > 3)
                            for st = 1 : length(locsF)
                                if st == 1
                                    if locsF(st)>wi+dm
                                        if (locsF(st+1)-locsF(st) > wi+dm)
                                            spkvr = horzcat(spkvr, Vr((locsF(st)-wi):(locsF(st)+wi)));
                                            spkvf = horzcat(spkvf, Vf((locsF(st)-wi):(locsF(st)+wi)));
                                            spkvs = horzcat(spkvs, Vs((locsF(st)-wi):(locsF(st)+wi)));
                                            spkvrrand = horzcat(spkvrrand,  Vr((locsRand(st)-wi):(locsRand(st)+wi)));
                                            spkvfrand = horzcat(spkvfrand,  Vf((locsRand(st)-wi):(locsRand(st)+wi)));
                                            spkvsrand = horzcat(spkvsrand,  Vs((locsRand(st)-wi):(locsRand(st)+wi)));
                                            spkdist = horzcat(spkdist, Vr(locsF(st)));
                                        end
                                    end
                                elseif st == length(locsF)
                                    if locsF(st)+wi+dm < length(Vr)
                                        if (locsF(st)-locsF(st-1) > wi+dm)
                                            spkvr = horzcat(spkvr, Vr((locsF(st)-wi):(locsF(st)+wi)));
                                            spkvf = horzcat(spkvf, Vf((locsF(st)-wi):(locsF(st)+wi)));
                                            spkvs = horzcat(spkvs, Vs((locsF(st)-wi):(locsF(st)+wi)));
                                            spkvrrand = horzcat(spkvrrand,  Vr((locsRand(st)-wi):(locsRand(st)+wi)));
                                            spkvfrand = horzcat(spkvfrand,  Vf((locsRand(st)-wi):(locsRand(st)+wi)));
                                            spkvsrand = horzcat(spkvsrand,  Vf((locsRand(st)-wi):(locsRand(st)+wi)));
                                            spkdist = horzcat(spkdist, Vr(locsF(st)));
                                        end
                                    end
                                else
                                    if (locsF(st)-locsF(st-1) > wi+dm && locsF(st+1)-locsF(st) > wi+dm)
                                        spkvr = horzcat(spkvr, Vr((locsF(st)-wi):(locsF(st)+wi)));
                                        spkvf = horzcat(spkvf, Vf((locsF(st)-wi):(locsF(st)+wi)));
                                        spkvs = horzcat(spkvs, Vs((locsF(st)-wi):(locsF(st)+wi)));
                                        spkvrrand = horzcat(spkvrrand,  Vr((locsRand(st)-wi):(locsRand(st)+wi)));
                                        spkvfrand = horzcat(spkvfrand,  Vf((locsRand(st)-wi):(locsRand(st)+wi)));
                                        spkvsrand = horzcat(spkvsrand,  Vs((locsRand(st)-wi):(locsRand(st)+wi)));
                                        spkdist = horzcat(spkdist, Vr(locsF(st)));
                                    end
                                end
                            end
                        end
                end
                SPKVR{n,j} = horzcat(SPKVR{n,j}, spkvr);
                SPKVF{n,j} = horzcat(SPKVF{n,j}, spkvf);
                SPKVS{n,j} = horzcat(SPKVS{n,j}, spkvs);
                SPKVRRand{n,j} = horzcat(SPKVRRand{n,j}, spkvrrand);
                SPKVFRand{n,j} = horzcat(SPKVFRand{n,j}, spkvfrand);
                SPKVSRand{n,j} = horzcat(SPKVSRand{n,j}, spkvsrand);
                SPKVRDist{n,j} = horzcat(SPKVRDist{n,j}, spkdist);
                SPKTIMES{n,j} = vertcat(SPKTIMES{n,j}, spikeTimes);
                VF{n,j} = vertcat(VF{n,j},Vf);
                VR{n,j} = vertcat(VR{n,j},Vr);
                VS{n,j} = vertcat(VS{n,j},Vs);
                ACST{n,j} = vertcat(ACST{n,j}, acSt);
            end
        end
        TOTALTIME(n,:) = totaltime;
    end
end
dtt.SPKVR = SPKVR;
dtt.SPKVF = SPKVF;
dtt.SPKVS = SPKVS;
dtt.SPKVRRand = SPKVRRand;
dtt.SPKVFRand = SPKVFRand;
dtt.SPKVSRand = SPKVSRand;
dtt.SPKVRDist = SPKVRDist;
dtt.SPKTIMES = SPKTIMES;
dtt.TOTALTIME = TOTALTIME;
dtt.VF = VF;
dtt.VR = VR;
dtt.VS = VS;
dtt.acSt = ACST;
pcav.EXP = cell(length(pTypes),1);
pcav.NSP = cell(length(pTypes),1);
pcav.PC1 = cell(length(pTypes),1);
pcav.SCO = cell(length(flies),length(pTypes));
pcav.Dist = cell(length(pTypes),1);
pcav.VrCents = -10:50:1500;
for j = 1 : length(pTypes)
    for nn = 1 : length(flies)
        
        [coeff,score,~,~,explained,~] = pca(SPKVR{nn,j}');
        if(size(coeff,2)>0)
            pcav.NSP{j} = vertcat(pcav.NSP{j}, size(SPKVR{nn,j},2));
            pcav.PC1{j} = horzcat(pcav.PC1{j}, coeff(:,1));
            pcav.SCO{nn,j} = horzcat(pcav.SCO{nn,j}, score(:,1));
            expl = zeros(6,1);
            expl(1:min(6,length(explained))) = explained(1:min(6,length(explained)));
            pcav.EXP{j} = horzcat(pcav.EXP{j}, expl);
            pcav.Dist{j} = vertcat(pcav.Dist{j}, hist(abs(SPKVRDist{nn,j}), ...
                pcav.VrCents)/sum(hist(abs(SPKVRDist{nn,j}), pcav.VrCents)));
        else
%             pcav.PC1{j} = horzcat(pcav.PC1{j}, zeros(41,1));
% %             pcav.SCO{nn,j} = horzcat(pcav.SCO{nn,j}, []);
%             expl = zeros(6,1);
%             expl(1:min(6,length(explained))) = explained(1:min(6,length(explained)));
%             pcav.EXP{j} = horzcat(pcav.EXP{j}, expl);
%             pcav.Dist{j} = vertcat(pcav.Dist{j}, hist(abs(SPKVRDist{nn,j}), ...
%                 pcav.VrCents)/sum(hist(abs(SPKVRDist{nn,j}), pcav.VrCents)));
        end
    end
end

end