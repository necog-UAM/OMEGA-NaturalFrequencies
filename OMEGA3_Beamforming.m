function OMEGA3_Beamforming (sub, ses, dpath)

% sub: subect code, e.g., '0001'
% ses: session, e.g., '0001'
% dpath: folder to store processed data

% -------------------------------------------------------- %
% 3.1. MRI normalization (output: mri_norm)
% 3.2. Head model
% 3.3. Forward model (output: source_forward_10mm)
% 3.4. Source reconstruction (output: source_inverse_10mm)
% -------------------------------------------------------- %

%% 3.1. MRI normalization
% mri normalized (mrin) and transformation matrix (normtrans)

cd([dpath '\sub-' sub '\ses-' ses])
load mri_coreg

cfg            = [];
cfg.nonlinear  = 'no';
mrin           = ft_volumenormalise(cfg, mri);      % do you want to change the anatomical labels for the axes [Y, n]? Y (r,a,s,i)
normtrans      = mrin.cfg.final;

save mri_norm mrin normtrans

%% 3.2. Head model
% semi-realistic singleshell head model based on the implementation from Guido Nolte

cfg             = [];
segment         = ft_volumesegment(cfg,mri);        % extract brain surface
segment.anatomy = mri.anatomy;

% Check that the segmentation is coregistered with mri
% figure
% cfg = [];
% cfg.interactive = 'yes';
% ft_sourceplot(cfg,mri);      % only mri
% cfg.funparameter = 'gray';
% ft_sourceplot(cfg,segment);  % segmented gray matter on top

cfg        = [];
cfg.method = 'singleshell';
vol        = ft_prepare_headmodel(cfg, segment);    % construct semi-realistic singleshell head model

%% 3.3. Forward model
% output: source_forward_10mm (contains mri, vol, grad, and the normalized grid and leadfields)  

cd([dpath '\sub-' sub '\ses-' ses])
load dataclean
grad  = dataclean.grad;

% Load normalized template grid (10mm)
load standard_sourcemodel3d10mm
grid = sourcemodel;

% Load normalized mri and head model
load mri_norm

% Adapt the normalized grid to each individual's brain space
posmni    = grid.pos;
pos       = ft_warp_apply(inv(normtrans), grid.pos*10, 'homogenous')/10;
grid.pos  = pos;
grid.unit = 'cm';

% Convert grad, vol and grid to common units (mm)
grad = ft_convert_units(grad, vol.unit);
grid = ft_convert_units(grid, vol.unit);

% Select only voxels within cortical mask (e.g. cerebellum is excluded)
% and corrected (inside the cortical surface projected with ft_sourceplot)
% created with select_corticalvox_aal
load ('D:\OMEGA-NaturalFrequencies\scripts\correccion_vox_inside_10mm.mat')
grid.inside = inside;

% Compute leadfields for each grid's voxel
cfg             = [];
cfg.grid        = grid;
cfg.grad        = grad;
cfg.vol         = vol;
cfg.channel     = {'MEG'};
cfg.normalize   = 'yes';
cfg.reducerank  = 2;
grid2           = ft_prepare_leadfield(cfg);

% Check that grad, vol and grid are correct (only for the first subject)
if strcmp(sub,'0001')
    figure
    plot3 (grad.chanpos(:,1), grad.chanpos(:,2), grad.chanpos(:,3), '.','MarkerEdgeColor',[0.8 0 0],'MarkerSize',25), hold on
    plot3 (vol.bnd.pos(:,1), vol.bnd.pos(:,2), vol.bnd.pos(:,3), '.','MarkerEdgeColor',[0 0 0.8]), hold on
    plot3 (grid2.pos(grid2.inside,1), grid2.pos(grid2.inside,2), grid2.pos(grid2.inside,3), '+k')
end

% Save grad, vol, grid and mri in source_forward structure to be used later
source_forward      = [];
source_forward.vol  = vol;
source_forward.mri  = mrin;
source_forward.grad = grad;
source_forward.grid = grid2;

save source_forward_10mm source_forward

%% 3.4. Computation of beamforming weights
% output: source_inverse_10mm (contains beamforming weights in source.avg.filter) 

cfg            = [];
cfg.covariance = 'yes';
datacov        = ft_timelockanalysis(cfg, dataclean);       % covariance matrix

% Compute spatial filters (in source.avg.filter)
cfg                   = [];
cfg.method            = 'lcmv';
cfg.grad              = source_forward.grad;
cfg.headmodel         = source_forward.vol;
cfg.grid              = source_forward.grid;
cfg.lcmv.fixedori     = 'yes';
cfg.lcmv.normalize    = 'yes';
cfg.lcmv.projectnoise = 'yes';
cfg.lcmv.keepfilter   = 'yes';          % important: save filters to use them later
cfg.lcmv.lambda       = '10%';          % the higher the smoother
cfg.lcmv.reducerank   = 2;
source                = ft_sourceanalysis(cfg, datacov);

load standard_sourcemodel3d10mm
source.avg.ori = {};
source.avg.mom = {};
source.avg.noisecov = {};
source.pos     = sourcemodel.pos;            % standard grid positions
source.inside  = grid.inside;

cd([dpath '\sub-' subs{s} '\ses-' sess{s}])
save source_inverse_10mm source

