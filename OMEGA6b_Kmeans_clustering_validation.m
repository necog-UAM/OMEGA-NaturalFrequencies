function OMEGA6b_Kmeans_clustering_validation 

% outpath: output folder
% replicate sample: 'tot' (N=128)
% cfg: specification of Nk (number of clusters to be tested), Nrep (number of repetitions of cross-validation procedure)

% -------------------------------------------------------- %
% 6b.1. Validation of the number of cluster by means of cross-validation
% 6b.2. Figures with results
% -------------------------------------------------------- %

%% 6b.1. Validation of the number of cluster by means of cross-validation

cd([outpath  'Nk25_10mm_' rep])
fil = sprintf('load kmeans_10mm_powsp_%s',rep);   
eval(fil)

Nobs  = size(powsptot,1);
Nks   = cfg.Nks;
Nrep  = cfg.Nrep;
Nperm = 100;                 % number of test datasets

D_perm = NaN(length(Nks),Nperm,Nrep);

for rep = 1:Nrep
    ct=1;
    for Nk = Nks
        disp(['Repetition ' num2str(rep) '/' num2str(Nrep) ' - Nclusters ' num2str(ct) '/' num2str(length(Nks))])
        r = randperm(Nobs);
        r_train = r(1:Nobs/2);
        r_rest  = r(Nobs/2+1:Nobs);

        powsptot_train = powsptot(r_train,:);
        [idx,C,sumd,D] = kmeans(powsptot_train,Nk,'Distance','cosine','Replicates',1,'MaxIter',200);

        Ntest = (Nobs/2) / Nperm;

        for i = 1:Nperm
            r_test = r_rest(Ntest*(i-1)+1:Ntest*i);
            powsptot_test = powsptot(r_test,:);
            [D_test,idx_test] = pdist2(C,powsptot_test,'cosine','Smallest',1);
            D_perm(ct,i,rep)  = sum(D_test);
        end
        ct=ct+1;
    end
end

save D_perm D_perm


%% 6b.2. Figures with results (Supplementary Figure 1A)

Nksx = Nks(1:end-1)+diff(Nks)/2;
col  = get(gca,'colororder');

% sumD
figure('Color', [1 1 1])
set(gcf,'Position',[0 0 600 250])
for rep = 1:Nrep
    plot(Nks,squeeze(mean(D_perm(:,:,rep),2)),'LineWidth',5,'Color',col(Nrep+1-rep,:)), hold on
end
set (gca,'XLim',[0 125],'YLim',[1.5e+4 3.5e+4],'FontName','Calibri','FontSize',30,'XTick',[25:25:100],'YTick',[20:20],'GridLineStyle','--','XGrid','on','GridAlpha',0.8)
box(gca,'off');
print('-dtiff','-r300','Nk_cross-validation_sumD.tiff');

% sumD slope
figure('Color', [1 1 1])
set(gcf,'Position',[0 0 600 250])
for rep = 1:Nrep
    plot(Nksx,squeeze(diff(mean(D_perm(:,:,rep),2))),'LineWidth',5,'Color',col(Nrep+1-rep,:)), hold on
end
set (gca,'XLim',[0 125],'FontName','Calibri','FontSize',30,'XTick',[25:25:100],'YTick',[20:20],'GridLineStyle','--','XGrid','on','GridAlpha',0.8)
box(gca,'off');
print('-dtiff','-r300','Nk_cross-validation_diffsumD.tiff');

