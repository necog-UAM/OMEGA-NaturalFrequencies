function OMEGA_SupplFig3_Single_subjects

cd([outpath  'Nk25_10mm_tot'])
eval(sprintf('load kmeans_10mm_powsp_tot'))   
eval(sprintf('load kmeans_10mm_Nk25_tot'))   
clear powsptot D

Nk = 25;
Nsub = 128;
badk = 21:25;

f   = 0.55:0.05:4.6;
foi = exp(f);

ff=[];
for k=1:Nk
    sp = C(k,:);
    f=(round(foi(find(sp==max(sp)))*10)/10);
    ff(k)=f;
end
[~,idf]=sort(ff);

ff = NaN(Nk,2);
for k = 1:Nk
    sp = C(idf(k),:);
    [pks,locs] = findpeaks(sp,'MinPeakHeight',0.1);
    
    if length(pks) == 1                               % unimodal spectrum
        ff(idf(k),1) = round(foi(locs(1))*10)/10;
        ff(idf(k),2) = round(foi(locs(1))*10)/10;     % repeat 1st peak
    elseif length(pks) == 2                           % bimodal spectrum
        ff(idf(k),1) = round(foi(locs(1))*10)/10;
        ff(idf(k),2) = round(foi(locs(2))*10)/10;
    end
end

propkz = (propk-mean(propk,2))./std(propk,0,2);       % z-normalize
propkz(idf(badk),:,:) = -10;                          % noisy clusters (any peak > 30Hz)
    
freqvoxs = [];
for s = 1:Nsub
    propkzm = propkz(:,:,s);
    for v = 1:Nvox
        fc = find(propkzm(:,v) == max(propkzm(:,v)));       % peak freq of the cluster with max Z-value
        if length(fc) > 1              % if more than one value, assign mode over neighbouring voxels
            fcneigh(1) = find(propkzm(:,v-1)==max(propkzm(:,v-1)));
            fcneigh(2) = find(propkzm(:,v-2)==max(propkzm(:,v-2)));
            fcneigh(3) = find(propkzm(:,v+1)==max(propkzm(:,v+1)));
            fcneigh(4) = find(propkzm(:,v+2)==max(propkzm(:,v+2)));
            fc = mode(fcneigh);
        end
        freqvoxs(v,s,:) = ff(fc,:);
    end
end

freqvoxs2 = freqvoxs(:,:,1);
for s = 1:Nsub
    for v = 1:Nvox
        if freqvoxs(v,s,1) ~= freqvoxs(v,s,2)
            freqvoxs2(v,s) = NaN ;
        end
    end
end

for s = 1:Nsub
    for v = 4:Nvox-3
        if isnan(freqvoxs2(v,s))
            ct=1;
            while isnan(freqvoxs2(v,s))
                fcneigh1(ct) = freqvoxs2(v-ct,s);
                fcneigh2(ct) = freqvoxs2(v+ct,s);
                ct = ct+1;
                fcneigh1(ct) = freqvoxs2(v-ct,s);
                fcneigh2(ct) = freqvoxs2(v+ct,s);
                ct = ct+1;
                fcneigh1(ct) = freqvoxs2(v-ct,s);
                fcneigh2(ct) = freqvoxs2(v+ct,s);
                freqvoxs2(v,s) =  mode([fcneigh1 fcneigh2]);
            end
        end
    end
end

load source_forward_10mm
load source_inverse_10mm
voxel_inside = find(source.inside==1);

cd(figpath)

for i=1:2
    figure('WindowState','maximized','Color',[1 1 1]);
    ct=1;
    for s = (i-1)*15+1:i*15
        source2 = source;
        source2.avg.pow(voxel_inside) = log(freqvoxs2(:,s));
        cfg=[];
        cfg.parameter  = 'avg.pow';
        cfg.downsample = 2;
        cfg.interpmethod = 'nearest';
        source_interp = ft_sourceinterpolate (cfg, source2, source_forward.mri);

        cfg               = [];
        cfg.method        = 'surface';
        cfg.funparameter  = 'pow';
        cfg.maskparameter = cfg.funparameter;
        cfg.funcolorlim   = [0.7 3.4];      % exponential scale   2-30 Hz
        cfg.funcolormap   = 'jet_omega_mod';
        cfg.projmethod    = 'nearest';
        cfg.opacity       = 0.8;
        cfg.camlight      = 'no';
        cfg.colorbar      = 'no';
        cfg.surffile     = 'surface_pial_right.mat';
        cfg.surfinflated  = 'surface_inflated_right_caret_white.mat';
        subplot(5,6,ct), ft_sourceplot(cfg,source_interp), view([90 0]),  camlight('right')
        ct=ct+1;
        subplot(5,6,ct), ft_sourceplot(cfg,source_interp), view([-90 0]), camlight('right')
        
        title(['sub' subs{s}])
        ct=ct+1;
    end
    print('-dtiff','-r300',['Single_subject_NFmaps' num2str(i) '_title.tiff']);
end

