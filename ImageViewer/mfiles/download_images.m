function download_images
% Download-images is an easy tool to download and create images for the
% later analysis

%load root_directory

params=parameter('Image downloader');

params=add(params,'directoryname','Select source directory name','./DCE');
params=add(params,'string','lesion name','lesion1');
params=add(params,'int','first slice','110');
params=add(params,'int','last slice','112');
params=add(params,'int','total slices/volume','224');

params=add(params,'button','Download','download_intern_detailed(params)');
params = add(params,'button','Quit','close');
params=parametergui(params);


