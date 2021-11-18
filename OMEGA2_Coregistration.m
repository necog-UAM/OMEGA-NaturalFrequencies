function OMEGA2_Coregistration (sub, ses, rawpath, dpath)

% sub: subect code, e.g., '0001'
% ses: session, e.g., '0001'
% rawpath: folder with raw data
% dpath: folder to store processed data

% -------------------------------------------------------- %
% 2.1. Coregistration of MEG-MRI spaces (output: mri_coreg)
% -------------------------------------------------------- %

%% 2.1. Coregistration of MEG-MRI spaces
% mri.transform is the transformation matrix to go from mri space to sensor space

cd([dpath '\sub-' sub])
mri = ft_read_mri('defaced_t1.nii');

cd([rawpath 'sub-' sub '\meg\meg2\omega\OMEGA_SCANS\sub-' sub '\ses-' ses '\meg\'])
dataresting = findfile('resting');
cd(dataresting)                                  

hsfile    = findfile('.pos');
headshape = ft_read_headshape(hsfile);    
[mri,scp] = omega_coreg([], mri, headshape);    % mark fiducials: lpa (l), rpa (r) and nasion (n), then quit (q) 

cd([dpath '\sub-' sub '\ses-' ses])
save mri_coreg mri scp

