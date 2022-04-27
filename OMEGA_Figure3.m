function OMEGA_Figure3 (outpath, rep, cfg)

% outpath: output folder
% replicate sample: 'tot' (N=128), 'rep1' (1st replicate sample N/2), 'rep2' (2nd replicate sample N/2)
% cfg: specification of Nk (number of clusters), Nvox (number of voxels), and Nsub (number of subjects)
% cfg.natfreq_fig : 1 to plot the brain map of natural frequencies in the cortical surface
% cfg.natfreq_nii : 1 to save nifti files

% -------------------------------------------------------- %
% F3.1. Plot the brain map of natural frequencies in the cortical surface
% F3.2. Save nifti files containing the brain map of natural frequencies
% -------------------------------------------------------- %

%% F3.1. Plot the brain map of natural frequencies in the cortical surface

Nk    = cfg.Nk;
Nsub  = cfg.Nsub;
Nvox  = cfg.Nvox;
cfg2  = cfg;

cd(outpath)
load source_forward_10mm
load source_inverse_10mm
voxel_inside = find(source.inside==1);
    
if strcmp(rep,'tot') || strcmp(rep,'rep1') || strcmp(rep,'rep2')
    cd([outpath 'Nk' num2str(Nk) '_10mm_' rep])
    load natfreq 
else
    cd([outpath 'Nk' num2str(Nk) '_10mm_' tot])
    load natfreq 
end

cd([outpath 'Nk' num2str(Nk) '_10mm_' rep])

if cfg2.natfreq_fig == 1
    source2 = source;
    if strcmp(rep,'tot') || strcmp(rep,'rep1') || strcmp(rep,'rep2')
        source2.avg.pow(voxel_inside) = log(natfmu);
    elseif strcmp(rep,'low95%CI') 
        source2.avg.pow(voxel_inside) = log(natfci(:,1));
    elseif strcmp(rep,'upp95%CI') 
        source2.avg.pow(voxel_inside) = log(natfci(:,2));
    end
    
    cfg=[];
    cfg.parameter  = 'avg.pow';
    cfg.downsample = 2;
    cfg.interpmethod = 'nearest';
    source_interp = ft_sourceinterpolate (cfg, source2, source_forward.mri);
        
    figure('WindowState','maximized','Color',[1 1 1]);

    cfg               = [];
    cfg.method        = 'surface';
    cfg.funparameter  = 'pow';
    cfg.maskparameter = cfg.funparameter;
    cfg.funcolorlim   = [0.7 3.4];      % exponential scale   2-30 Hz
    cfg.funcolormap   = 'jet_omega_mod';   
    cfg.projmethod    = 'nearest';
    cfg.opacity       = 0.8;
    cfg.camlight      = 'no';
    cfg.colorbar      = 'yes';
    cfg.surffile     = 'surface_pial_left.mat';
    cfg.surfinflated  = 'surface_inflated_left_caret_white.mat';
    subplot(2,2,1), ft_sourceplot(cfg,source_interp), view([-90 0]), camlight('left')
    subplot(2,2,3), ft_sourceplot(cfg,source_interp), view([90 0]),  camlight('left')
    
    cfg.surffile     = 'surface_pial_right.mat';
    cfg.surfinflated  = 'surface_inflated_right_caret_white.mat';
    subplot(2,2,2), ft_sourceplot(cfg,source_interp), view([90 0]),  camlight('right')
    subplot(2,2,4), ft_sourceplot(cfg,source_interp), view([-90 0]), camlight('right')
    
    print('-dtiff','-r300',['Nk' num2str(Nk) '_Natural_Frequencies.tiff']);
end



%% F3.2. Save nifti files containing the brain map of natural frequencies and 95% CI

if cfg2.natfreq_nii == 1
    
    cd([outpath 'Nk' num2str(Nk) '_10mm_' rep])
    format shortG             
    
    source2 = source;
    if strcmp(rep,'tot') || strcmp(rep,'rep1') || strcmp(rep,'rep2')
        source2.avg.pow(voxel_inside) = round(natfmu*10)/10;             % round to one decimal place in nifti files
    elseif strcmp(rep,'low95%CI')
        source2.avg.pow(voxel_inside) = round(natfci(:,1)*10)/10;
    elseif strcmp(rep,'upp95%CI')
        source2.avg.pow(voxel_inside) = round(natfci(:,2)*10)/10;
    end
    
    cfg=[];
    cfg.parameter  = 'avg.pow';
    cfg.downsample = 2;
    cfg.interpmethod = 'nearest';
    source_interp = ft_sourceinterpolate (cfg, source2, source_forward.mri);
    
    cfg = [];
    cfg.filetype  = 'nifti';
    cfg.parameter = 'pow';
    if strcmp(rep,'tot') || strcmp(rep,'rep1') || strcmp(rep,'rep2')
        cfg.filename  = 'NaturalFreq.nii';
    elseif strcmp(rep,'low95%CI')
        cfg.filename  = 'NaturalFreq_Low95%CI.nii';
    elseif strcmp(rep,'upp95%CI')
        cfg.filename  = 'NaturalFreq_High95%CI.nii';
    end
    ft_sourcewrite(cfg, source_interp)
    
    natf = spm_vol(cfg.filename);
    aal = spm_vol([outpath 'raal_cortical_mask.nii']);
    natf.fname = cfg.filename;
    dat = natf.private.dat(:,:,:).*aal.private.dat(:,:,:);          % masked with cortical mask
    spm_write_vol(natf,dat)
    
end



