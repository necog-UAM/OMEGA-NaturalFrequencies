function OMEGA7_Naturalfreq (outpath, rep, cfg)

% outpath: output folder
% replicate sample: 'tot' (N=128), 'rep1' (1st replicate sample N/2), 'rep2' (2nd replicate sample N/2)
% cfg: specification of Nk (number of clusters), Nvox (number of voxels), Nsub (number of subjects), and Nboot (number of bootstrap repetitions
% cfg.bootci: 1 to compute natural frequencies for each voxel and bootstrap CIs

% -------------------------------------------------------- %
% 7.1. Compute natural frequency of each voxel (output: natfreq)
% 7.2. Compute natural frequency for each AAL region (output: natfreq_aal)
% -------------------------------------------------------- %

%% 7.1. Compute natural frequency of each voxel (output: natfreq)

Nk    = cfg.Nk;
Nsub  = cfg.Nsub;
Nvox  = cfg.Nvox;
Nboot = cfg.Nboot;
Nperm = cfg.Nperm;
badk  = cfg.badk;
cfg2  = cfg;

if strcmp(rep,'rep1') || strcmp(rep,'rep2')
    Nsub = cfg.Nsub/2;
end

cd(outpath)
load source_forward_10mm
load source_inverse_10mm
voxel_inside = find(source.inside==1);

if cfg2.bootci == 1                          
    cd([outpath  'Nk' num2str(Nk) '_10mm_' rep])
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
        
    ff = NaN(Nk,2);
    for k = 1:Nk
        sp = C(idf(k),:);
        [pks,locs] = findpeaks(sp,'MinPeakHeight',0.1);
        
        if length(pks) == 1                                % unimodal spectrum
            ff(idf(k),1) = round(foi(locs(1))*10)/10;
            ff(idf(k),2) = round(foi(locs(1))*10)/10;      % if unimodal spectrum, repeat 1st peak
        elseif length(pks) == 2                            % bimodal spectrum
            ff(idf(k),1) = round(foi(locs(1))*10)/10;
            ff(idf(k),2) = round(foi(locs(2))*10)/10;
        end
    end
    
    propkz = (propk-mean(propk,2))./std(propk,0,2);        % z-normalize
    propkz(idf(badk),:,:) = -10;                           % noisy clusters (any peak > 30Hz)

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
            freqvoxs(v,s,:) = ff(fc);
        end
    end
    
    if cfg.bimod ==1        % check if each voxel came from unimodal or bimodal spectra
        for v=1:Nvox
            modalf(v) = 100*sum((freqvoxs(v,:,1)-freqvoxs(v,:,2))~=0)./Nsub;
        end    

        % Figure 4
        source2 = source;
        source2.avg = rmfield(source2.avg,'mom');
        source2.avg.pow(voxel_inside) = modalf;

        cfg = [];
        cfg.parameter  = 'avg.pow';
        cfg.downsample = 2;
        cfg.interpmethod  =  'linear';
        source_interp  = ft_sourceinterpolate (cfg, source2, source_forward.mri);

        figure('Color',[1 1 1]);
        set(gcf,'Position',[100 100 1100 900])

        cfg               = [];
        cfg.method        = 'surface';
        cfg.funparameter  = 'pow';
        cfg.maskparameter = cfg.funparameter;
        cfg.funcolorlim   = [0  35];
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

        print('-dtiff','-r300',['Figure4_bimodality.tiff']);
    end
    
    
    f1 = 0;
    f2 = 4.6;
    d  = 0.2;
    foi1 = exp(f1:d:f2);
    foi2 = exp(f1:d/10:f2);
    for i = 1:length(foi2)
        c0 = foi2(i);
        rectp(i,:) = rectangularPulse(c0-c0/4 ,c0+c0/4 ,foi2);      % to isolate individual peaks later      
    end
    
    p = parpool(4);                  % parallel pool
    p.IdleTimeout = inf;
    
    natf.A       = NaN(Nvox,Nboot);
    natf.mu      = NaN(Nvox,Nboot);
    natf.sigma   = NaN(Nvox,Nboot);
    natf.rsq     = NaN(Nvox,Nboot);
    natf.randsub = NaN(Nboot,Nsub);
    
    % bootstrap confidence interval; based on bootci -> ci = bootci(Nboot,{@bootfun,freqvoxs},'type','per')
    rng('shuffle')
    for b = 1:Nboot
        disp(['Bootstrapping ' num2str(b) '/' num2str(Nboot)])
        randsub = randi(Nsub,Nsub);
        randsub = randsub(1,:);                       % sample with replacement Nsub
        
        for v = 1:Nvox
            fv = freqvoxs(v,randsub,:);
            [counts,~] = histcounts(fv(:),foi1);
            centers = exp(log(foi1) + d/2);
            counts  = interp1(centers(1:end-1),counts,foi2,'pchip');
            centers = foi2;
                
            [pks,locs] = findpeaks(counts);
            [~,ii] = sort(pks,'descend');
            if length(ii) >= 3
                c0s = centers(locs(ii(1:3)));         % fit a maximum of 3 peaks to speed up
            else
                c0s = centers(locs);
            end
                                   
            % gaussian fit of all candidate peak frequencies
            [A,mu,sigma,rsq] = gausfitc0(c0s,counts,centers,rectp,foi2);
               
            i = find(rsq == max(rsq));          % identify the one with the highest goodness of fit
            
            natf.A(v,b)       = A(i);
            natf.mu(v,b)      = mu(i);
            natf.sigma(v,b)   = sigma(i);
            natf.rsq(v,b)     = rsq(i);
            natf.randsub(b,:) = randsub;
            
            gausf = A(i) * exp(-(centers-mu(i)).^2 /(2*sigma(i).^2));
            plot(centers,counts,'k'), hold on
            plot(centers,squeeze(gausf),'r')
            set(gca,'XLim',[0 35],'YLim',[0 80])
            hold off
            pause
        end
    end
    
    natfmu      = median(natf.mu,2);               
    natfci(:,1) = prctile(natf.mu,2.5,2);          % 95% CI interval  
    natfci(:,2) = prctile(natf.mu,97.5,2);
    
    cd([outpath  'Nk' num2str(Nk) '_10mm_' rep])
    save natfreq natf natfmu natfci  
    
    
    if cfg.rand_ci == 1                            % statistical test of 95% CIs
        
        natfmu = NaN(Nperm,1);
        natfci = NaN(Nperm,2);
        
        for pp = 1:Nperm
            disp(['Permutation ' num2str(pp) '/' num2str(Nperm)])
            natf.mu      = NaN(Nboot,1);
            
            for b = 1:Nboot
                randsub = randi(Nsub,Nsub);
                randsub = randsub(1,:);                       % sample with replacement Nsub
                
                randvox  = randperm(Nvox);
                freqvoxs = freqvoxs(randvox,:,:);
                
                for v = 1
                    fv = freqvoxs(v,randsub,:);
                    [counts,~] = histcounts(fv(:),foi1);
                    centers = exp(log(foi1) + d/2);
                    counts  = interp1(centers(1:end-1),counts,foi2,'pchip');
                    centers = foi2;
                    
                    [pks,locs] = findpeaks(counts);
                    [~,ii] = sort(pks,'descend');
                    if length(ii) >= 3
                        c0s = centers(locs(ii(1:3)));         % fit a maximum of 3 peaks to speed up
                    else
                        c0s = centers(locs);
                    end
                    
                    % gaussian fit of all candidate peak frequencies
                    [A,mu,sigma,rsq] = gausfitc0(c0s,counts,centers,rectp,foi2);
                    
                    i = find(rsq == max(rsq));          % identify the one with the highest goodness of fit
                    
                    natf.mu(b)      = mu(i);
                    
                    gausf = A(i) * exp(-(centers-mu(i)).^2 /(2*sigma(i).^2));
                end
            end
            
            natfmu(pp) = median(natf.mu);
            natfci(pp,1) = prctile(natf.mu,2.5);          % 95% CI interval
            natfci(pp,2) = prctile(natf.mu,97.5);
        end
        
        cd([outpath  'Nk' num2str(Nk) '_10mm_' rep])
        fil = sprintf('save natfreq_randci natfmu natfci',Nk);
        eval(fil)
        
    end
    
    delete(p)
end


%% 7.2. Compute natural frequency for each AAL region (output: natfreq_aal)

cd([outpath  'Nk' num2str(Nk) '_10mm_' rep])
load natfreq 
load([outpath 'aal_voxel_label_10mm.mat'])
Nroi = length(aal_label_reduc);

for roi = 1:Nroi
    inds = find(label_inside_aal_reduc==roi);
    voxs = voxel_inside_aal(inds);
    muvoxs = natf.mu(voxs,:);
    
    Nk = 4;
    [idx,C,sumd,D] = kmeans(muvoxs,Nk,'Replicates',5,'MaxIter',200);      % k-means 4 clusters max
    nvoxs = [sum(idx==1) sum(idx==2) sum(idx==3) sum(idx==4) ];           % number of voxels in each cluster                       
    Nk = Nk - sum(nvoxs==1);                                              % eliminate clusters with only 1 member and repeat clustering                                           
    [idx,C,sumd,D] = kmeans(muvoxs,Nk,'Replicates',5,'MaxIter',200);      % k-means adapting number of clusters 
    nvoxs = [sum(idx==1) sum(idx==2) sum(idx==3) sum(idx==4) ];           % number of voxels in each cluster                       

    dominclust = find(nvoxs == max(nvoxs));                               % dominant cluster
    centroidvox = find(D(:,dominclust)==min(D(:,dominclust)));            % voxel closer to centroid
    represvox(roi) = voxs(centroidvox(1));                                % representative voxel for this AAL region
%      histogram(C'),title(aal_label_reduc{roi}),pause                    % visualize pattern of natural frequencies for both clusters
end

natfmu_aal      = median(natf.mu(represvox,:),2);
natfci_aal(:,1) = prctile(natf.mu(represvox,:),2.5,2);       % 95% CI interval
natfci_aal(:,2) = prctile(natf.mu(represvox,:),97.5,2);
% [natfmu_aal natfci_aal]

cd([outpath 'Nk25_10mm_' rep])
save natfreq_aal natfmu_aal natfci_aal represvox aal_label_reduc

