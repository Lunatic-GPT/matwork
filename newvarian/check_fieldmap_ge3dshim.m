function check_fieldmap_ge3dshim

a=read_fdf('B0.f.fdf');
write_afni(a*62,'fieldmap');
