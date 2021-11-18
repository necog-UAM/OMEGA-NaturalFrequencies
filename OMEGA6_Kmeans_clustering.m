function OMEGA6_Kmeans_clustering (outpath, rep, cfg)

% outpath: output folder
% replicate sample: 'tot' (N=128), 'rep1' (1st replicate sample N/2), 'rep2' (2nd replicate sample N/2)
% cfg: specification of Nk (number of clusters), Nvox (number of voxels), Ntr (number of time segments), and Nsub (number of subjects)

% -------------------------------------------------------- %
% 6.1. Computation of K-means clustering
% 6.2. Proportion of power spectra (out of 150) categorized in each cluster
% -------------------------------------------------------- %

%% 6.1. Computation of K-means clustering

% Read matrix with the power spectra of all subjects, time segments and voxels
cd([outpath  'Nk25_10mm_' rep])
fil = sprintf('load kmeans_10mm_powsp_%s',rep);   
eval(fil)

% K-means clustering
Nk = cfg.Nk;
[idx,C,sumd,D] = kmeans(powsptot,Nk,'Distance','cosine','Display','iter','Replicates',5,'MaxIter',200);


%% 6.2. Proportion of power spectra (out of 150) categorized in each cluster
% This computation is done for each subject and voxel (for each
% subject and voxel, the sum across clusters is equal to 1)

Nvox = cfg.Nvox;
Ntr  = cfg.Ntr;
Nsub = cfg.Nsub;

propk = NaN(Nk,Nvox,Nsub);

for k = 1:Nk
    disp(['Cluster ' num2str(k) '/' num2str(Nk)])
    ctk = find(idx==k);
    ctksub = ksub(ctk)';
    ctkvox = kvox(ctk)';
    ctksubfilt = ctksub==[1:Nsub];
    ctkvoxfilt = ctkvox==[1:Nvox];
    for s=1:Nsub
        propk(k,:,s)=sum(ctksubfilt(:,s) & ctkvoxfilt)./Ntr;
    end
end

cd([outpath  'Nk25_10mm_' rep])
fil = sprintf('save kmeans_10mm_Nk%d_%s idx C sumd D Nvox propk -v7.3',Nk,rep)
eval(fil)

