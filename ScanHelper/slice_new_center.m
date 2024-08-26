function slice_new_center(center,an,shift)
% input:
% center: old center [y,z]; dicom convention
% an: angle in degrees: should be T2C
% shift: how much to shift in mm; + for upward shift

if ~exist('shift','var')
    shift=15;
end

vec=[-sin(an*pi/180),cos(an*pi/180)];
res=center+shift*vec;

fprintf('new center:');
disp(res);

