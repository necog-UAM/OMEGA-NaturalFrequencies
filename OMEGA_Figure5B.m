function OMEGA_Figure5B (outpath, figpath, rep, cfg)

% outpath: output folder
% figpath: folder to save figure
% replicate sample: 'tot' (N=128), 'rep1' (1st replicate sample N/2), 'rep2' (2nd replicate sample N/2)
% cfg: specification of Nk (number of clusters), Nvox (number of voxels), and Nsub (number of subjects)

% -------------------------------------------------------- %
% F5B.1. Proportion of time in each frequency band and AAL region
% -------------------------------------------------------- %

%% F5B.1. Proportion of time in each frequency band and AAL region

Nk   = cfg.Nk;
Nsub = cfg.Nsub;
Nvox = cfg.Nvox;

cd([outpath  'Nk25_10mm_' rep])
fil = sprintf('load kmeans_10mm_Nk%d_%s',Nk,rep);
eval(fil)

cd([outpath 'Nk25_10mm_' rep])
load natfreq_aal
load([outpath 'aal_voxel_label_10mm.mat'])
aal  = spm_vol([outpath 'raal_cortical_mask.nii']);
aal  = aal.private.dat(:,:,:);
Nroi = length(aal_label_reduc);

f   = 0.55:0.05:4.6;
foi = exp(f);
ff  = [];
for k = 1:Nk
    sp = C(k,:);
    f  = round(foi(find(sp==max(sp)))*10)/10;
    ff(k) = f;
end
[~,idf] = sort(ff);

Nk2   = sum(ff < 30);   % select only clusters with peak centroids <30Hz (>30Hz are probably artefacts) 
foi   = foi(1:findbin(foi,35));

fbands{1} = [1 2 3];                  % delta
fbands{2} = [4 5 6];                  % theta
fbands{3} = [7 8  10 11  13 14];      % alpha
fbands{4} = [9 12];                   % mu
fbands{5} = [15 16];                  % beta low
fbands{6} = [17 18 19];               % beta high 

Nfb = length(fbands);

for roi = 1:Nroi
    p = propkm(:,represvox(roi));
    for fb = 1:Nfb
        pfbands(roi,fb) = sum(p(idf(fbands{fb})));       
    end
end

figure('Color',[1 1 1]);
set(gcf,'Position',[100 0 500 1000])
imagesc(pfbands)
colormap(flipud(gray));        % change colormap to gray (higher values are black)ff

textStrings = num2str(pfbands(:), '%0.1f');       
textStrings = strtrim(cellstr(textStrings));  
[x,y] = meshgrid(1:Nfb,1:Nroi); 
hStrings = text(x(:),y(:),textStrings(:),'HorizontalAlignment','center','FontSize',8);
midValue = mean(get(gca, 'CLim'));                
textColors = repmat(pfbands(:) > midValue, 1, 3);       % choose white or black for the text color of the strings depending on background colour
set(hStrings, {'Color'}, num2cell(textColors, 2));  
set(gca, 'XTick', 1:Nfb, 'XTickLabel', {'delta', 'theta', 'alpha', 'mu', 'beta low','beta high'}, ...
         'YTick', 1:Nroi, 'YTickLabel', aal_label_reduc, 'TickLength', [0 0], 'FontSize',11);
           
cd(figpath)
print('-dtiff','-r300',['Proptime_aal_fbands.tiff']);

