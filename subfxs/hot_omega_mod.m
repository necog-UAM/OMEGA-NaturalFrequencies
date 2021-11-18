function J = hot_omega_mod
%JET    Variant of HSV
%   JET(M) returns an M-by-3 matrix containing the jet colormap, a variant
%   of HSV(M). The colors begin with dark blue, range through shades of
%   blue, cyan, green, yellow and red, and end with dark red. JET, by
%   itself, is the same length as the current figure's colormap. If no
%   figure exists, MATLAB uses the length of the default colormap.
%
%   See also PARULA, HSV, HOT, PINK, FLAG, COLORMAP, RGBPLOT.

%   Copyright 1984-2015 The MathWorks, Inc.

load('hotxmod.mat');
J = cmap;
