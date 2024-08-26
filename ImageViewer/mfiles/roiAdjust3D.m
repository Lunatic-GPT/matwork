function roiAdjust3D



analysis=parameter('roiAdjust3D');

analysis=add(analysis,'directoryname','Select data directory name',fullfile('./','Data'));
analysis=add(analysis,'int','first slice','1');
analysis=add(analysis,'int','total number of slices','3');
analysis=add(analysis,'button','Run','roiAdjust3D_callback(params)');
analysis=add(analysis,'button','Quit','close');


parametergui(analysis);
