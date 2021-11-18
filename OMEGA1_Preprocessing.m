function OMEGA1_Preprocessing (sub, ses, rawpath, dpath)

% sub: subect code, e.g., '0001'
% ses: session, e.g., '0001'
% rawpath: folder with raw data
% dpath: folder to store processed data

% -------------------------------------------------------- %
% 1.1. Reading data
% 1.2. Preprocessing (output: data)
% 1.3. Artifact correction with ICA (output: dataclean, comp)
% -------------------------------------------------------- %

%% 1.1. Reading data

cd([rawpath 'sub-' sub '/meg/meg2/omega/OMEGA_SCANS/sub-' sub '/ses-' ses '/meg/'])

dataresting = findfile('resting');
cd(dataresting)                                  

datafile    = findfile('.meg4');
headerfile  = findfile('.res4');

cfg            = [];
cfg.datafile   = datafile;
cfg.headerfile = headerfile;
cfg.continuous = 'yes';
data           = ft_preprocessing(cfg);

%% 1.2. Preprocessing 

cfg         = [];
cfg.detrend = 'yes';
data        = ft_preprocessing(cfg,data);

cfg            = [];
cfg.refchannel = 'MEGREF';
data           = ft_denoise_pca(cfg,data);

cfg            = [];
cfg.hpfilter   = 'yes';
cfg.hpfreq     = 0.05;
cfg.hpfiltord  = 3;
data           = ft_preprocessing(cfg,data);
data.trial{1}  = ft_preproc_dftfilter(data.trial{1},data.fsample,[60 120 180],'dftreplace','neighbour');    % Mewett et al., 2004

cfg            = [];
cfg.resamplefs = 512;
cfg.detrend    = 'yes';
cfg.demean     = 'yes';
data           = ft_resampledata(cfg,data);

cfg         = [];
cfg.toilim  = [10 data.time{1}(end)-10];
data        = ft_redefinetrial(cfg,data);

cd([dpath '/sub-' sub '/ses-' ses])
save data data

%% 1.3. Artifact correction with ICA 

cfg              = [];
cfg.numcomponent = 40;                          % keep 40 principal components               
cfg.method       = 'runica';
comp             = ft_componentanalysis(cfg ,data);         

cd([dpath '/sub-' sub '/ses-' ses])
save comp comp

% Detect badicas and identify badsegments (many ICs contaminated during the same time segment)
cfg        = [];
cfg.layout = 'CTF275.lay';
cfg.scale  = 4e-11; 
badicas    = detect_badicas_omega (cfg, data, comp);      % automatically save badicas: mouse-click reject; space-bar keep component
% badicas = [];
% save badicas badicas                                    % manually save badicas 

% Remove badicas from signal and save dataclean
cfg           = [];
cfg.component = badicas;
dataclean     = ft_rejectcomponent (cfg, comp, data);

cd([dpath '/sub-' sub '/ses-' ses])
save dataclean dataclean

