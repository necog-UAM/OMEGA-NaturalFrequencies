function [mri,scp] = omega_coreg(cfg, mri, shape)

% Contributors: Gavin Paterson, Jan-Mathijs Schoffelen, Almudena Capilla & Joachim Gross
% University of Glasgow

if ~isfield(cfg, 'ne'),          cfg.ne        = 1;          end
if ~isfield(cfg, 'ng'),          cfg.ng        = 1;          end
if ~isfield(cfg, 'threshold'),   cfg.threshold = [0.05 0.5]; end
if ~isfield(cfg, 'verify'),      cfg.verify      = 'yes';    end

%load the data if necessary
if isstr(mri),   mri   = read_fcdc_mri(mri);    end;
if isstr(shape), shape = read_headshape(shape); end;

%get shape in correct units
shape  = ft_convert_units(shape, 'mm');
nshape = size(shape.pos,1);

%first do an approximate manual coregistration
tmpcfg        = [];
tmpcfg.method = 'interactive';
try, mri.transformorig = mri.transform; end
mri.transform = eye(4);
mri           = ft_volumerealign(tmpcfg, mri);

%extract scalp surface
tmpcfg = [];
tmpcfg.thr.low  = cfg.threshold(1);
tmpcfg.thr.high = cfg.threshold(2);
tmpcfg.ne       = cfg.ne;
tmpcfg.ng       = cfg.ng;
scpvox          = scalp3(mri.anatomy, tmpcfg);
%[scpvox,scptri]  = scalp3b(mri.anatomy);

%get into head coordinates according to the crude coregistration
scp     = ft_warp_apply(mri.transform, scpvox);
scporig = scp;
scp     = scp(1:5:end,:);

%do icp FIXME try to understand the last input argument
%FIXME2 there's no reason to tease apart the translation and rotation see checkdata
[R, t, corr, D, data2] = icp2(scp', shape.pos', 20);
M             = [R t;0 0 0 1];
mri.transformorig1 = mri.transform;
mri.transform = inv(M)*mri.transform;
scp           = ft_warp_apply(mri.transform, scpvox);

if strcmp(cfg.verify, 'yes'),
    figure;hold on;
    plot3(scp(1:5:end,1),scp(1:5:end,2),scp(1:5:end,3),'.','markersize',4);
    plot3(shape.pos(:,1),shape.pos(:,2),shape.pos(:,3),'r.');
    try,plot3(shape.fid.pos(1:3,1),shape.fid.pos(1:3,2),shape.fid.pos(1:3,3),'go','markersize',6,'linewidth',4); end
    try,
        fiduc2 = mri.cfg.fiducial;
        fiduc2 = [fiduc2.nas;fiduc2.lpa;fiduc2.rpa];
        fiduc2 = ft_warp_apply(mri.transformorig1,fiduc2);
        plot3(fiduc2(:,1),fiduc2(:,2),fiduc2(:,3),'ko','markersize',6,'linewidth',4);
    end
    view([0 0]);axis vis3d;axis off
end



