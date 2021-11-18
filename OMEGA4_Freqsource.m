function OMEGA4_Freqsource (sub, ses, dpath)

% sub: subect code, e.g., '0001'
% ses: session, e.g., '0001'
% dpath: folder to store processed data

% -------------------------------------------------------- %
% 4.1. Reconstruction of source-level time-series
% 4.2. Frequency analysis parameters
% 4.3. Frequency analysis computation (output: freq_allvox_10mm)
% -------------------------------------------------------- %

%% 4.1. Reconstruction of source-level time-series

% Read necessary files
cd([dpath '\sub-' sub '\ses-' ses])
load dataclean
load source_forward_10mm
load source_inverse_10mm

% Reconstruct source-space data
time         = dataclean.time{1};
voxel_inside = find(source.inside==1);
Nvox         = length(voxel_inside);
datasource   = zeros(Nvox,length(time));
for i = 1:Nvox
    disp([num2str(i) ' / ' num2str(length(voxel_inside))])
    datasource(i,:) = source.avg.filter{voxel_inside(i)} * dataclean.trial{1};
end
datasource = datasource./repmat(std(datasource,0,2),[1,length(time)]);        % baseline correction to account for the centre of the head bias

% Discard artifactual time segments 
cd([dpath '\sub-' sub '\ses-' ses])
if exist('badsegments.mat') == 2
    load badsegments
    for b = 1:length(badsegments)
        t1 = findbin(time,badsegments{b}(1));
        t2 = findbin(time,badsegments{b}(2));
        time(t1:t2) = [];
        datasource(:,t1:t2) = [];
    end
end

% Use only 5-minute recordings for all participants
if time(end) > 290      % >5 min
    t1 = findbin(time,290);
    t2 = findbin(time,time(end));
    time(t1:t2) = [];
    datasource(:,t1:t2) = [];
end

% Organize source-reconstructed data in a new fieldtrip structure
dataclean.trial{1}   = single(datasource);
dataclean.time{1}    = time;
dataclean.label      = {};
dataclean.sampleinfo = [1 size(dataclean.trial{1},2)];
for i = 1:Nvox
    dataclean.label{i} = ['V' num2str(i)];
end
clear datasource

% If there were not badsegments, data is organized in 1 single trial
% If there were badsegments (discontinuities), data is organized in several trials
if exist('badsegments.mat')==2
    dtime      = diff(time);
    [pks,locs] = findpeaks(dtime);
    pkstime    = time(locs);
    if length(pkstime) > 0
        trl = [];
        trl(1,1) = 1;
        for t = 1:length(pks)
            tt = findbin(time, pkstime(t));
            trl(t,2) = tt - 1;
            trl(t+1,1) = tt + 1;
        end
        trl(t+1,2) = length(time);
        trl(:,3)   = 0;
        
        cfg     = [];
        cfg.trl = trl;
        dataclean = ft_redefinetrial(cfg,dataclean);     
    end
end


%% 4.2. Frequency analysis parameters: foi, toi, t_ftimwin

f         = 0.55:0.05:4.6;
foi       = exp(f);                % logarithmically spaced frequencies
t_ftimwin = 5./foi;              % time-window adapted to each frequency (5 cycles)
t_fstep   = 0.5;                   % sliding time-window (500 ms steps)

bd   = t_ftimwin(1).*2;       % remove borders = 2*time-window
bdpt = bd.*dataclean.fsample;

Ntr  = 0;
toi2 = {};
ct   = 1;
valid_tr = [];
for tr = 1:length(dataclean.trial)
    toi = bd:t_fstep:dataclean.time{tr}(end)-bd;
    if length(toi) > 0
        toi2{ct} = toi;
        Ntr = Ntr+length(toi);
        valid_tr(ct) = tr;
        ct = ct+1;
    end
end


%% 4.3. Frequency analysis computation

Nvox  = length(dataclean.label);
Nf    = length(foi);
powsp = single(NaN(Nvox,Nf,Ntr));

tr2 = 0;
for tr = 1:length(valid_tr)
    tr1 = tr2+1;
    tr2 = tr1+length(toi2{tr})-1;
    
    cfg            = [];
    cfg.trials     = valid_tr(tr);
    cfg.method     = 'mtmconvol';
    cfg.taper      = 'hanning';
    cfg.output     = 'pow';
    cfg.foi        = foi;
    cfg.toi        = toi2{tr};
    cfg.t_ftimwin  = t_ftimwin;
    cfg.pad        = 'nextpow2';
    cfg.keeptrials = 'yes';
    freq           = ft_freqanalysis(cfg, dataclean);
    
    powsp(:,:,tr1:tr2)= single(squeeze(freq.powspctrm));
end

cd([dpath '\sub-' sub '\ses-' ses])
save freq_allvox_10mm powsp foi

