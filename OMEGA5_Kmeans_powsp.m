function datasegments = OMEGA5_Kmeans_powsp (subs, sess, dpath, outpath, rep)

% subs: all subjects
% sess: corresponding sessions
% dpath: folder with partipants data
% outpath: output folder
% replicate sample: 'tot' (N=128), 'rep1' (1st replicate sample N/2), 'rep2' (2nd replicate sample N/2)
% datasegments: number of valid time segments per participant

% -------------------------------------------------------- %
% 5.1. Select sample: N or N/2
% 5.2. Preparation of power spectra matrix for subsequent K-means clustering (output: kmeans_10mm_powsp_tot/rep1/rep2)
% -------------------------------------------------------- %

%% 5.1. Select sample: N or N/2

cd(outpath)
load replicates_goodsubj
replicates = {'rep1','rep2','tot'};   

if find(strcmp(replicates,rep)) == 1       
    goodsubj = goodsubj1;      
elseif find(strcmp(replicates,rep)) == 2
    goodsubj = goodsubj2;      
elseif find(strcmp(replicates,rep)) == 3  
    goodsubj = 1:length(subs);      
end 

Nsub = length(goodsubj);
 
%% 5.2. Preparation of power spectra matrix for subsequent K-means clustering

Ntr      = 150;
powsptot = [];
ksub     = [];
kvox     = [];
ktrial   = [];
rng('shuffle')
datasegments = [];

for s = 1:Nsub
    disp(['Subject ' num2str(s) '/' num2str(Nsub)])
    cd([dpath '\sub-' subs{s} '\ses-' sess{s}])
    load freq_allvox_10mm 
    ct = 1;
    valid_tr = [];
    for tr = 1:size(powsp,3)
        if sum(sum(isnan(powsp(:,:,tr)))) == 0
           valid_tr(ct) = tr; 
           ct = ct+1;
        end
    end
    datasegments(s) = length(valid_tr);          % number of valid time segments per subject

    rndtr   = randperm(length(valid_tr),Ntr);
    select_tr = valid_tr(rndtr);
    powsp2  = powsp(:,:,select_tr);
    ct = 1;
    powsp3  = single(zeros(size(powsp2,1)*size(powsp2,3),size(powsp2,2)));
    ksub2   = [];
    kvox2   = [];
    ktrial2 = [];
    for i = 1:size(powsp2,1)
        for tr = 1:size(powsp2,3)
            powsp3(ct,:) = powsp2(i,:,tr);
            ksub2(ct)    = s;
            kvox2(ct)    = i;
            ktrial2(ct)  = select_tr(tr);
            ct = ct+1;
        end
    end
    powsptot = [powsptot; powsp3];         % concatenation of power spectra
    ksub     = [ksub ksub2];               % keep track of subject, voxel and trial of each power spectrum
    kvox     = [kvox kvox2];
    ktrial   = [ktrial ktrial2];  
end

clear freq powsp

bl       = sum(powsptot,2);                                % save kmeans_10mm_baseline bl
powsptot = powsptot./repmat(bl,[1 size(powsptot,2)]);      % compute relative power to correct the center of the head bias

cd([outpath  'Nk25_10mm_' rep])
fil = sprintf('save kmeans_10mm_powsp_%s powsptot ksub kvox ktrial -v7.3',rep);
eval (fil)


