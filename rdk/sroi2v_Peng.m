p =parameter('convert surface rois to volume');
p = add(p,'string','subject','BQ');
p = add(p,'filename','aligned local anatomy','BQ_SurfVol_ns_Alnd_Exp+orig');

p = add(p,'directoryname','roi directory',pwd);
p = add(p,'string','left surface rois','V1.lh.1D.roi,V2dx.lh.1D.roi,V2v.lh.1D.roi,V3dx.lh.1D.roi,V3v.lh.1D.roi,V4.lh.1D.roi,MTLoc.lh.1D.roi');
p = add(p,'string','right surface rois','V1.rh.1D.roi,V2d.rh.1D.roi,V2v.rh.1D.roi,V3d.rh.1D.roi,V3v.rh.1D.roi,V4.rh.1D.roi,MTLoc.rh.1D.roi');

p = add(p,'directoryname','suma','../../SUMA');

p = add(p,'filename','brain mask','brainmask+orig');

sroi2v_gui_callback(p);
