function OMEGA_Figure2 (outpath, figpath, rep, cfg)

% outpath: output folder
% figpath: folder to save localization of generators in cortical surface (F1.4.)
% replicate sample: 'tot' (N=128), 'rep1' (1st replicate sample N/2), 'rep2' (2nd replicate sample N/2)
% cfg: specification of Nk (number of clusters), Nvox (number of voxels), and Nsub (number of subjects)
% cfg.fig_powsp: 1 to plot centroid power spectra 
% cfg.stats: 1 to compute stats;
% cfg.fig_brainpatt: 1 to plot brain distribution associated to each cluster
% cfg.locgenerators_vox: list of peak voxels generating different brain rhythms
% cfg.locgenerators_fig: 1 to plot the localization of generators in cortical surface
% cfg.fig_gamma: 1 to plot (artefactual) generators of gamma-band activity

% -------------------------------------------------------- %
% F2.1. Plot centroid power spectra for each cluster (1-30Hz)
% F2.2. Brain pattern of each cluster: stats
% F2.3. Brain pattern of each cluster: figure (masked with T-stat)
% F2.4. Localize generators in cortical surface
% F2.5. Supplementary Figure 1: plot high-frequency (>30 Hz) generators in slices
% -------------------------------------------------------- %

%% F2.1. Plot centroid power spectra for each cluster (1-30Hz)

Nk   = cfg.Nk;
Nsub = cfg.Nsub;
Nvox = cfg.Nvox;
cfg2 = cfg;

cd([outpath  'Nk25_10mm_' rep])
fil = sprintf('load kmeans_10mm_Nk%d_%s',Nk,rep);
eval(fil)

f   = 0.55:0.05:4.6;
foi = exp(f);
ff  = [];
for k = 1:Nk
    sp = C(k,:);
    f  = round(foi(find(sp==max(sp)))*10)/10;
    ff(k) = f;
end
[~,idf]=sort(ff);

Nk2   = sum(ff < 30);   % select only clusters with peak centroids <30Hz (>30Hz are most probably artefacts) 
foi   = foi(1:findbin(foi,35));
c     = jet_omega_mod;
fqlog = exp(0.7:0.044:3.5);

for k = 1:Nk2
    sp = C(idf(k),1:findbin(foi,35));
    [pks,locs] = findpeaks(sp,'MinPeakHeight',0.1);
    if length(locs) == 2     % if mu rhythm, then assign the frequency of the beta peak
        locs = locs(2);
    end
    f  = round(foi(locs)*10)/10;
    fc = findbin(fqlog,f);
    
    if cfg2.fig_powsp == 1
        figure('Color', [1 1 1])
        set(gcf,'Position',[100 100 550 350])
        
        plot(foi,sp,'Color',c(fc,:),'LineWidth',6)
        title(['k ' num2str(k),' - ',num2str(f),' Hz'])
        set (gca,'XLim',[1.5 35],'YLim',[0 .4],'FontName','Calibri','FontSize',30,'XTick',[10 20 30],'YTickLabel',{'','','','',''})
        box(gca,'off');
        print('-dtiff','-r300',['Nk' num2str(Nk) '_Power_Spectrum_k' num2str(k),' - ',num2str(f), 'Hz.tiff']);
    end
end


%% F2.2. Brain pattern of each cluster: stats

cd(outpath)
load source_forward_10mm
load source_inverse_10mm
voxel_inside = find(source.inside==1);

propkz = (propk-mean(propk,2))./std(propk,0,2);        % Z-normalize
propkzm = mean(propkz,3);

if cfg2.stats == 1
    for k = 1:Nk
        disp(['Cluster ' num2str(k) '/' num2str(Nk)])
        d_obs  = squeeze(propkz(k,:,:))';
        d_obsm = repmat(mean(d_obs,2),[1,Nvox]);
        [h,p,ci,stats] = ttest(d_obs,d_obsm);
        t_obs(k,:) = stats.tstat;
        for it = 1:1000
            for s = 1:Nsub
                vox_shuf = randperm(Nvox);
                dx = squeeze(propkz(k,vox_shuf,s))';
                d_shuf(s,:)= dx;
            end
            d_shufm = repmat(mean(d_shuf,2),[1,Nvox]);
            [h,p,ci,stats] = ttest(d_shuf,d_shufm);
            t_max = max(stats.tstat);
            t_shuf(it) = t_max;
        end
        t_thr05(k)  = prctile(t_shuf,95);
        t_thr01(k)  = prctile(t_shuf,99);
        t_thr001(k) = prctile(t_shuf,99.9);
    end
    
    cd([outpath  'Nk25_10mm_' rep])
    save kmeans_stats_DEF idf t_obs t_thr05 t_thr01 t_thr001
end


%% F2.3. Brain pattern of each cluster: figure (masked with T-stat)

if cfg2.fig_brainpatt == 1
    cd([outpath  'Nk25_10mm_' rep])
    load kmeans_stats_DEF
    
    for k = 1:Nk
        t_sig = t_obs(idf(k),:);
        t_sig(t_sig < t_thr05(idf(k)))  = 0;
        t_sig(t_sig >= t_thr05(idf(k))) = 1;
        
        sp = C(idf(k),1:findbin(foi,30));
        [pks,locs] = findpeaks(sp,'MinPeakHeight',0.1);
        if length(locs) == 2          % if mu rhythm, then assign the frequency of the beta peak
            locs = locs(2);
        end
        f = round(foi(locs)*10)/10;
        
        source2 = source;
        source2.avg.pow(voxel_inside) = propkzm(idf(k),:).*t_sig;        % Z-value masked by stat
        
        temp = source2.avg.pow(voxel_inside);                 % check minimum Z-value significant
        temp(temp ==0) = NaN;
        zthr(k)=nanmin(temp);
        
        cfg = [];
        cfg.parameter  = 'avg.pow';
        cfg.downsample = 2;
        cfg.interpmethod  =  'linear';
        source_interp  = ft_sourceinterpolate (cfg, source2, source_forward.mri);
        
        % Write MRIcroN image
        cfg = [];
        cfg.filetype  = 'nifti';
        cfg.parameter = 'pow';
        cfg.filename  = ['Nk' num2str(Nk) '_k' num2str(k) '_' num2str(round(f)) 'Hz'];
        ft_sourcewrite(cfg, source_interp)
        
        datx = spm_vol([cfg.filename '.nii']);
        aal  = spm_vol([outpath 'raal_cortical_mask.nii']);
        datx.fname = [cfg.filename '.nii'];
        dat = datx.private.dat(:,:,:).*aal.private.dat(:,:,:);
        spm_write_vol(datx,dat)
        
        % Plot and save results in cortical surface
        figure('Color',[1 1 1]);
        set(gcf,'Position',[100 100 1100 900])
        
        cfg               = [];
        cfg.method        = 'surface';
        cfg.funparameter  = 'pow';
        cfg.maskparameter = cfg.funparameter;
        cfg.funcolorlim   = [0  1.2];
        cfg.funcolormap   = 'hot_omega_mod';
        cfg.projmethod    = 'nearest';
        cfg.opacitymap    = 'rampup';
        cfg.camlight      = 'no';
        cfg.colorbar      = 'yes';
        cfg.surffile      = 'surface_pial_left.mat';
        cfg.surfinflated  = 'surface_inflated_left_caret_white.mat';
        subplot(2,2,1), ft_sourceplot(cfg,source_interp), view([-90 0]), camlight('left')
        subplot(2,2,3), ft_sourceplot(cfg,source_interp), view([90 0]),  camlight('left')
        
        cfg.surffile      = 'surface_pial_right.mat';
        cfg.surfinflated  = 'surface_inflated_right_caret_white.mat';
        subplot(2,2,2), ft_sourceplot(cfg,source_interp), view([90 0]),  camlight('right')
        subplot(2,2,4), ft_sourceplot(cfg,source_interp), view([-90 0]), camlight('right')
        
        print('-dtiff','-r300',['Nk' num2str(Nk) '_k' num2str(k) '_' num2str(round(f)) 'Hz.tiff']);
    end
end


%% F2.4. Localize generators in cortical surface

if cfg2.locgenerators_fig == 1
    cd(figpath)
    vx = cfg2.locgenerators_vox;
    for i = 1:length(vx)
        v(i) = find(voxel_inside==vx(i));
    end
    
    for i = 1:length(v)
        figure('Color',[1 1 1]);
        set(gcf,'Position',[100 100 1100 900])
        
        source2 = source;
        source2.avg.pow(voxel_inside) = 0;        
        source2.avg.pow(voxel_inside(v(i))) = 1;    
        
        cfg=[];
        cfg.parameter  = 'avg.pow';
        cfg.downsample = 2;
        cfg.interpmethod  = 'nearest';
        source_interp = ft_sourceinterpolate (cfg, source2, source_forward.mri);
        
        cfg               = [];
        cfg.method        = 'surface';
        cfg.funparameter  = 'pow';
        cfg.maskparameter = cfg.funparameter;
        cfg.funcolormap   = 'hot_omega_mod';   % parula
        cfg.projmethod    = 'nearest';
        cfg.opacitymap    = 'rampup';
        cfg.camlight      = 'no';
        cfg.colorbar      = 'yes';
        cfg.surffile      = 'surface_pial_left.mat';
        cfg.surfinflated  = 'surface_inflated_left_caret_white.mat';
        subplot(2,2,1), ft_sourceplot(cfg,source_interp), view([-90 0]), camlight('left')
        subplot(2,2,3), ft_sourceplot(cfg,source_interp), view([90 0]),  camlight('left')
        
        cfg.surffile      = 'surface_pial_right.mat';
        cfg.surfinflated  = 'surface_inflated_right_caret_white.mat';
        subplot(2,2,2), ft_sourceplot(cfg,source_interp), view([90 0]),  camlight('right')
        subplot(2,2,4), ft_sourceplot(cfg,source_interp), view([-90 0]), camlight('right')
        
        print('-dtiff','-r300',['AAL_voxel' num2str(vx(i)) '.tiff']);
        close(gcf)
    end
end


%% F2.5. Supplementary Figure 1: plot high-frequency (>30 Hz) generators in slices

if cfg2.fig_gamma == 1
    
    Nk   = cfg.Nk;
    Nsub = cfg.Nsub;
    Nvox = cfg.Nvox;
    cfg2 = cfg;
    
    cd([outpath  'Nk25_10mm_' rep])
    fil = sprintf('load kmeans_10mm_Nk%d_%s',Nk,rep);
    eval(fil)
    
    f   = 0.55:0.05:4.6;
    foi = exp(f);
    ff  = [];
    for k = 1:Nk
        sp = C(k,:);
        f  = round(foi(find(sp==max(sp)))*10)/10;
        ff(k) = f;
    end
    [~,idf]=sort(ff);
    
    cd(outpath)
    load source_forward_10mm
    load source_inverse_10mm
    voxel_inside = find(source.inside==1);
    
    propkz = (propk-mean(propk,2))./std(propk,0,2);        % Z-normalize
    propkzm = mean(propkz,3);
    
    cd([outpath  'Nk25_10mm_' rep])
    load kmeans_stats_DEF
    
    for k = 20:Nk
        t_sig = t_obs(idf(k),:);
        t_sig(t_sig < t_thr05(idf(k)))  = 0;
        t_sig(t_sig >= t_thr05(idf(k))) = 1;
        
        sp = [C(idf(k),:) 0];
        [pks,locs] = findpeaks(sp,'MinPeakHeight',0.1);
        if length(locs) == 2          % if mu rhythm, then assign the frequency of the beta peak
            locs = locs(2);
        end
        f = round(foi(locs)*10)/10;
        
        source2 = source;
        source2.avg.pow(voxel_inside) = propkzm(idf(k),:).*t_sig;        % Z-value masked by stat
        
        temp = source2.avg.pow(voxel_inside);                 % check minimum Z-value significant
        temp(temp ==0) = NaN;
        zthr(k)=nanmin(temp);
        
        cfg = [];
        cfg.parameter  = 'avg.pow';
        cfg.downsample = 2;
        cfg.interpmethod  =  'linear';
        source_interp  = ft_sourceinterpolate (cfg, source2, source_forward.mri);
        
        figure('Color',[1 1 1]);
        set(gcf,'Position',[100 100 1100 900])
        cfg               = [];
        cfg.method        = 'slice';
        cfg.nslices       = 10;
        cfg.funparameter  = 'pow';
        cfg.funcolorlim   = [0 2.5];
        cfg.maskparameter = cfg.funparameter;
        cfg.funcolormap   = 'hot';
        cfg.projmethod    = 'nearest';
        % cfg.interactive   = 'yes';
        cfg.inputcoord    = 'mni';
        ft_sourceplot(cfg,source_interp)
        
        print('-dtiff','-r300',['SupplFig1_Nk' num2str(Nk) '_k' num2str(k) '_' num2str(round(f)) 'Hz.tiff']);
        
    end
end

