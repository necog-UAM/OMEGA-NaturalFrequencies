function [nf, nf_low, nf_upp ,aal] = NaturalFreq (mnicoord, dpath)

% mnicoord: MNI coordinates for a given voxel (1x3 vector)
% dpath: folder containing nifti files and templates:
% - 'NaturalFreq.nii'
% - 'NaturalFreq_Low95%CI.nii'
% - 'NaturalFreq_High95%CI.nii'
% - 'VoxelsInside_10mm.nii'
% - 'aal_voxel_label_10mm.mat'
% nf: natural frequency of the voxel
% nf_low: lower bound of the 95% CI 
% nf_upp: upper bound of the 95% CI 
% aal: AAL region corresponding to the MNI coordinate

if exist('dpath') == 0
    dpath = pwd;
end
cd(dpath)

d = spm_vol ('voxels_inside_10mm.nii');
coord = [];
for i = 1:3
    c0 = -d.mat(i,4)./d.mat(i,i);
    coord(i) = c0 + mnicoord (i)./d.mat(i,i);
end
vox = d.private.dat(coord(1),coord(2),coord(3));

% Find natural frequency and 95% CI
d      = spm_vol ('NaturalFreq.nii');
nf     = d.private.dat(coord(1),coord(2),coord(3));

d      = spm_vol ('NaturalFreq_Low95%CI.nii');
nf_low = d.private.dat(coord(1),coord(2),coord(3));

d      = spm_vol ('NaturalFreq_High95%CI.nii');
nf_upp = d.private.dat(coord(1),coord(2),coord(3));

% Find AAL region
load('aal_voxel_label_10mm.mat');
v = find(voxel_inside_aal == vox);

% Display results
if ~isempty(v) && nf > 0 
    aal = aal_label{label_inside_aal(v)};
    disp(['Natural Frequency at MNI coordinates (' num2str(mnicoord(1)) ',' num2str(mnicoord(2)) ',' num2str(mnicoord(3)) ...
    '): '  num2str(nf) ' Hz [' num2str(nf_low) ' - ' num2str(nf_upp) '] - ' aal])
else
    disp(['No cortical region found at MNI coordinates (' num2str(mnicoord(1)) ',' num2str(mnicoord(2)) ',' num2str(mnicoord(3)) ')'])
end


%% Create template for voxels_inside

% load source_forward_10mm
% load source_inverse_10mm
% voxel_inside = find(source.inside==1);
% 
% source.avg.pow(voxel_inside) = [1:length(voxel_inside)];
% cfg              = [];
% cfg.parameter    = 'pow';
% cfg.downsample   = 2;
% cfg.interpmethod = 'nearest';
% template         = ft_sourceinterpolate(cfg, source, source_forward.mri);
% 
% cfg = [];
% cfg.filetype  = 'nifti';
% cfg.parameter = 'pow';
% cfg.filename  = ['VoxelsInside_10mm.nii'];
% ft_sourcewrite(cfg, template)

