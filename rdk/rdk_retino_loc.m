sid = 'K004';
sv = sprintf('%s_SurfVol_ns_Alnd_Exp+orig',sid);
suma = '../../segmentation/SUMA';
cmask = '-a brainmask+orig -expr step(a)';

vol2surf('cc_gamma_coh2+orig',sid,cmask,suma,sv);
correlateAnalysis('ts_sort_coh3_surf.lh.niml.dset','ref_gamma.1D','cc_coh3.lh');
correlateAnalysis('ts_sort_coh3_surf.rh.niml.dset','ref_gamma.1D','cc_coh3.rh');

%vol2surf('rdkave_vr_dtr+orig',sid,cmask,suma,sv);
correlateAnalysis('rdkave_vr_dtr_surf.lh.niml.dset','ref_gamma.1D','cc.lh');
correlateAnalysis('rdkave_vr_dtr_surf.rh.niml.dset','ref_gamma.1D','cc.rh');

vol2surf('ringave_vr+orig',sid,cmask,suma,sv);
vol2surf('wedgesave_vr+orig',sid,cmask,suma,sv);
retino_wdgrng('wedgesave_vr_surf.lh.niml.dset','ringave_vr_surf.lh.niml.dset','retino.lh',16);
retino_wdgrng('wedgesave_vr_surf.rh.niml.dset','ringave_vr_surf.rh.niml.dset','retino.rh',16);