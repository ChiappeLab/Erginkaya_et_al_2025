function [VForwSeg, pTypes] = GetForwSegVData(path, params)
flies = dir(path);
flies = flies(3:end);
dt = load([path flies(1).name '\DataLowRes.mat'], 'Flies');
if contains(path,'Dark')
    temp = dt.Flies.Seq;
    for i = 1:length(temp)
        temp{i} = 'Dark';
    end
    dt.Flies.Seq = temp;
    clearvars temp
end
seq = dt.Flies.Seq;
pTypes = unique(seq);
VForwSeg = cell(length(pTypes),length(flies));
for n = 1 : length(flies)
    pathF = [path flies(n).name];
    pathi = [pathF '\DataLowRes.mat'];
    if exist(pathi, 'file') == 2
        dt = load(pathi);
        dt = dt.Flies;
        seq = dt.Seq;
        if contains(pathi,'Dark')
            temp = seq;
            for i = 1:length(temp)
                temp{i} = 'Dark';
            end
            seq = temp;
            clearvars temp
        end
        pTypes = unique(seq);
        dt = dt.Data;
        for k = 1 : length(dt)
            for l = 1 : length(pTypes)
                switch seq{k}
                    case pTypes{l}
                        vrb = dt{k}.Vr;
                        vsb = dt{k}.Vs;
                        vfb = dt{k}.Vf;
                        actst = dt{k}.actState;
                        [locs, ~, cmhSong, thr] = SortSpikes(vrb, actst, params);
                        [locsF, pksIF,pksEF] = CullSpikesBasedOnVf(vrb, locs, cmhSong, thr, vfb, params);
                        [locsF, pksIF, pksEF] = TemplateCompSpike(vrb,locsF,pksIF, pksEF, cmhSong, thr, params);
                        [FBouts] = GetFbouts(actst, vfb, locsF, pksIF,pksEF, params);
                       
                        cSac = length(VForwSeg{l,n});
                        for i = 1 : length(FBouts)
                            if (FBouts{i}(1)) > 0 && (FBouts{i}(end)) < length(dt{k}.FramesC2)
                                frms = (dt{k}.FramesC2(FBouts{i}(1)):dt{k}.FramesC2(FBouts{i}(end)))-dt{k}.FramesC2(1);
                                if min(frms) > 0
                                    VForwSeg{l,n}{i+cSac}.VrLR = vrb(FBouts{i});
                                    VForwSeg{l,n}{i+cSac}.VfLR = vfb(FBouts{i});
                                    VForwSeg{l,n}{i+cSac}.VsLR = vsb(FBouts{i});
                                end
                            end
                        end
                end
            end
        end
    end
    disp(['Fly ',num2str(n),'/',num2str(length(flies)),' Done'])
    disp(pathF)
end
end
