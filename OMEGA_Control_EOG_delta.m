function [R_blinks, P_blinks, R_saccs, P_saccs] = OMEGA_Control_EOG_delta (subs, sess, rawpath, dpath, outpath)    

% Control analysis to test the potential correlation between delta activity (z-values) in the 
% orbitofrontal cortex and both blinks and saccadic movements estimated from EOG recordings

%% Identification of blinks/saccades in EOG traces

Nsub   = length(subs);
blinks = NaN(Nsub,1);
saccs  = NaN(Nsub,1);

for s = 1:Nsub
    sub = subs{s};
    ses = sess{s};
    
    dataresting=[];
    try
        cd([rawpath 'sub-' sub '/meg/meg2/omega/OMEGA_SCANS/sub-' sub '/ses-' ses '/meg/'])
        dataresting = findfile('resting');
    end
    try
        cd([rawpath 'sub-' sub '/ses-' ses '/meg/'])
        dataresting = findfile('.ds');
    end
    cd(dataresting)

    datafile    = findfile('.meg4');
    headerfile  = findfile('.res4');

    cfg            = [];
    cfg.datafile   = datafile;
    cfg.headerfile = headerfile;
    cfg.continuous = 'yes';
    data           = ft_preprocessing(cfg);
  
    chveog1 = [];
    chveog2 = [];
        
    chheog1 = [];
    chheog2 = [];
    
    for i = 1:length(data.label)
        ch = strfind(data.label{i},'VEOG');
        if length(ch) == 1
            chveog1 = i;
        end
    end
    
    if length(chveog1) == 1
        chheog1 = match_str(data.label,'HEOG');
    elseif length(chveog1) == 0
        chveog2 = match_str(data.label,'EEG058');
        chheog2 = match_str(data.label,'EEG059');
    end
    
    if length(chveog1) == 1 & length(chheog1) == 1
        cfg  = [];
        cfg.channel = {'VEOG','HEOG'};
    elseif length(chveog1) == 1 & length(chheog1) == 0
        cfg  = [];
        cfg.channel = {'VEOG'};
    end
    
    if length(chveog2) == 1 & length(chheog2) == 1
        cfg  = [];
        cfg.channel = {'EEG058','EEG059'};
    elseif length(chveog2) == 1 & length(chheog2) == 0
        cfg  = [];
        cfg.channel = {'EEG058'};
    end
    
    data           = ft_preprocessing(cfg,data);
    
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

    time         = data.time{1};
    if exist('badsegments.mat') == 2
        load badsegments
        for b = 1:length(badsegments)
            t1 = findbin(time,badsegments{b}(1));
            t2 = findbin(time,badsegments{b}(2));
            time(t1:t2) = [];
            data.trial{1}(:,t1:t2) = NaN;
        end
    end
    
    if time(end) > 290      % >5 min
        t1 = findbin(time,290);
        t2 = findbin(time,time(end));
        time(t1:t2) = [];
        data.trial{1}(:,t1:t2) = NaN;
    end

    try, veog{s} = data.trial{1}(1,:); end
    try, heog{s} = data.trial{1}(2,:); end
    
    try, blinks(s) = nanstd(veog{s});  end
    try, saccs(s)  = nanstd(heog{s});  end
      
end

% noisy VEOG channels
badveog = [66];
% noisy HEOG channels
badheog = [16 117];

blinks2 = blinks;
saccs2  = saccs;

blinks2(badveog) = NaN;
saccs2(badheog)  = NaN;


%% Identification of the orbitofrontal voxel with maximal activity in delta clusters

cd(outpath)
load kmeans_10mm_Nk25_tot
load source_inverse_10mm
voxel_inside = find(source.inside==1);
vx = find(voxel_inside==2925);         % voxel in orbitofrontal cortex with max z-value

f   = 0.55:0.05:4.6;
foi = exp(f);

Nk = 25;
ff = [];
for k = 1:Nk
    sp = C(k,:);
    f  = round(foi(find(sp==max(sp)))*10)/10;
    ff(k) = f;
end
[~,idf] = sort(ff);
 
delta_ks = idf(1:2);                                     % delta clusters
propkz   = (propk-mean(propk,2))./std(propk,0,2);        % z-normalize

delta_z  = squeeze(mean(propkz(delta_ks,vx,:),1));

save eog_delta veog heog delta_z blinks2 saccs2

%% Correlation coefficient between delta_z and eye movements (blinks/saccades)

[R_blinks, P_blinks] = corrcoef(delta_z,blinks2,'rows','pairwise')
[R_saccs , P_saccs]  = corrcoef(delta_z,saccs2,'rows','pairwise')

