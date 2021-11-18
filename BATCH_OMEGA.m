%% ----------- General information ---------------- %

dpath = 'D:\OMEGA-NaturalFrequencies\data\';
rawpath = 'D:\OMEGA-NaturalFrequencies\raw\';       

sub = [1 2 3 4 5 6 7 8 9 11 12 14 15 16 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 37 39 40 41 42 44 45 46 47 48 49 50 51 52 55 56 57 58 59 60 61 62 63 64 65 67 68 69 70 71 72 73 74 75 76 77 78 79 80 84 85 87 88 89 90 91 92 94 95 96 97 98 99 101 102 103 104 105 106 134 145 146 148 149 150 151 152 154 155 156 157 158 159 160 161 165 166 167 168 169 170 171 175 176 177 179 181 184 185 195 197 200 207 208 210 212]';
ses = [1 1 1 1 1 1 1 2 1  2  1  2  1  1  1  2  3  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  3  1  2  1  2  1  2  1  1  2  1  1  1  1  1  1  1  1  2  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 ]';

subs = {} ; sess = {};
for s = 1:length(sub)
    if sub(s) < 10
        subs{s} = ['000' num2str(sub(s))];
    elseif sub(s) >= 10 & sub(s) < 100 
        subs{s} = ['00' num2str(sub(s))];
    else
        subs{s} = ['0' num2str(sub(s))];
    end
    sess{s} = ['000' num2str(ses(s))];
    
end

Nsub   = length(sub);

%% ----------- Reading, preprocessing and artifact correction ---------------- %
  
for s = 1:Nsub
    OMEGA1_Preprocessing (subs{s}, sess{s}, rawpath, dpath)     
end

%% ----------- Source reconstruction (computation of beamforming weights) ---------------- %
 
% MEG-MRI coregistration: mark fiducials (lpa, rpa and nasion; then quit) 
% and check result (if coregistration is not accurate, repeat procedure)
for s = 1:Nsub
    OMEGA2_Coregistration (subs{s}, sess{s}, rawpath, dpath)         
end

% Compute forward model (source_forward_10mm) and beamforming weights (source_inverse_10mm) (10mm-grid)
for s = 1:Nsub
    OMEGA3_Beamforming (subs{s}, sess{s}, dpath)           
end

%% ----------- Frequency analysis on reconstructed source-space time-series ---------------- %

for s = 1:Nsub
    OMEGA4_Freqsource (subs{s}, sess{s}, dpath)         
end

%% ----------- K-means clustering of source-space power spectra ---------------- %

% Preparation of power spectrum matrix as needed for K-means clustering: 36960000 (128 subjs x 150 time-segments x 1925 voxels)  X  82 freqs
% Repeat analysis for the complete sample ('tot', N=128), 'rep1', and 'rep2' (random subsamples of N/2)

replicates = {'rep1','rep2','tot'};   
outpath    = 'D:\OMEGA-NaturalFrequencies\results\';      
replicate_randsamples(outpath)              % to obtain 2 random subsamples of 64 participants

rep = 'tot';  OMEGA5_Kmeans_powsp (subs, sess, dpath, outpath, rep)         
rep = 'rep1'; OMEGA5_Kmeans_powsp (subs, sess, dpath, outpath, rep)         
rep = 'rep2'; OMEGA5_Kmeans_powsp (subs, sess, dpath, outpath, rep)         

% Computation of K-means clustering (it requires near 32 Gb of RAM for the complete sample)
% K-means parameters: 'Distance' = 'cosine', 'Replicates' = 5, 'MaxIter' = 200
% Repeat analysis for the whole sample ('tot', N=128), 'rep1', and 'rep2' (random subsamples of N/2)

cfg = [];
cfg.Nk   = 25;
cfg.Nvox = 1925;
cfg.Ntr  = 150;
cfg.Nsub = 128;
rep = 'tot';  OMEGA6_Kmeans_clustering (outpath, rep, cfg)         

cfg.Nsub = 128/2;
rep = 'rep1'; OMEGA6_Kmeans_clustering (outpath, rep, cfg)         
rep = 'rep2'; OMEGA6_Kmeans_clustering (outpath, rep, cfg)         


%% ---------------------- FIGURE 2 --------------------------- %
% Plots of centroid power spectra and brain distribution for each cluster (1-30Hz) 

replicates = {'rep1','rep2','tot'};   
outpath    = 'D:\OMEGA-NaturalFrequencies\results\';      
figpath    = 'D:\OMEGA-NaturalFrequencies\figures\AAL_regions';
rep = 'tot';  

cfg = [];
cfg.Nk   = 25;
cfg.Nvox = 1925;
cfg.Nsub = 128;
cfg.fig_powsp = 1;
cfg.stats     = 0;
cfg.fig_brainpatt = 0;
cfg.locgenerators_fig = 0;
cfg.locgenerators_vox = [2479 2472 2925 4816 4636 3578 3929 2066 1717 4637 4634 5053 4634 5126 5117 4261 4633 5070 5047 5863 5477 4459 4451 3722 4091];     % list of generators of brain rhythms
cfg.fig_gamma = 0;          % Supplementary Figure 1

OMEGA_Figure2 (outpath, figpath, rep, cfg) 

%% ----------- Computation of natural frequencies for each voxel and AAL region ---------------- %

replicates = {'rep1','rep2','tot'};   
outpath    = 'D:\OMEGA-NaturalFrequencies\results\';      

cfg        = [];
cfg.Nk     = 25;
cfg.Nvox   = 1925;
cfg.Nsub   = 128;
cfg.Nboot  = 500;       
cfg.bootci = 1;

rep = 'tot';  OMEGA7_Naturalfreq (outpath, rep, cfg) 
rep = 'rep1'; OMEGA7_Naturalfreq (outpath, rep, cfg) 
rep = 'rep2'; OMEGA7_Naturalfreq (outpath, rep, cfg) 


% Matlab function to look up natural frequencies for given MNI coordinates
rep = 'tot';
cd([outpath  'Nk25_10mm_' rep])
mnicoord = [2,-62,38];
[nf, nf_low, nf_upp ,aal] = NaturalFreq (mnicoord);

%% ---------------------- FIGURES 3-4 --------------------------- %
% Plots of brain maps of natural frequencies for the whole sample and replicate samples

replicates = {'rep1','rep2','tot','low95%CI','upp95%CI'};   
outpath    = 'D:\OMEGA-NaturalFrequencies\results\';      

cfg       = [];
cfg.Nk    = 25;
cfg.Nvox  = 1925;
cfg.Nsub  = 128;
cfg.natfreq_fig = 0;
cfg.natfreq_nii = 1;

rep = 'tot';      OMEGA_Figure3 (outpath, rep, cfg) 
rep = 'rep1';     OMEGA_Figure3 (outpath, rep, cfg) 
rep = 'rep2';     OMEGA_Figure3 (outpath, rep, cfg) 
rep = 'low95%CI'; OMEGA_Figure3 (outpath, rep, cfg) 
rep = 'upp95%CI'; OMEGA_Figure3 (outpath, rep, cfg) 


%% ----------- Correlation between replicate1 and replicate2 brain maps ---------------- %

cd([outpath 'Nk25_10mm_rep1'])
load natfreq
natfmu1 = natfmu;

cd([outpath 'Nk25_10mm_rep2'])
load natfreq
natfmu2 = natfmu;

[R,P] = corrcoef(natfmu1,natfmu2)


%% ---------------------- FIGURE 5 --------------------------- %
% Spectral modes: z-values for each AAL region and frequency band

replicates = {'rep1','rep2','tot'};   
outpath    = 'D:\OMEGA-NaturalFrequencies\results\';      
figpath    = 'D:\OMEGA-NaturalFrequencies\figures\';
rep = 'tot';

cfg = [];
cfg.Nk   = 25;
cfg.Nvox = 1925;
cfg.Nsub = 128;

OMEGA_Figure5A (outpath, figpath, rep, cfg)
OMEGA_Figure5B (outpath, figpath, rep, cfg)
