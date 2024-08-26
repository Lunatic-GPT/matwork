function auto_shift_DCE_images
% auto_shift_DCE_images 
% automatically align the images slice by slice with 2dImReg function in
% afni.
% 5/5/2009, Wrote it. XZ 

%load root_directory
shifty=parameter('DCE image shifting (x/y shifts)');

shifty=add(shifty,'directoryname','Select image directory name',fullfile('./','DCE_temp'));
shifty=add(shifty,'directoryname','Select data directory name',fullfile('./','Data'));
%shifty=setcallback(shifty,'Select image directory name', 'save shift_gui');


%shifty=add(shifty,'string','Enter the subject ID','N84432');
%shifty=setcallback(shifty,'Enter the subject ID','save shift_gui');
shifty=add(shifty,'int','First slice','1');

shifty=add(shifty,'int','Total number of slices','3');
%shifty=setcallback(shifty,'Total number of slices', 'save shift_gui');
shifty=add(shifty,'bool','Perform motion correction',true);
%shifty=add(shifty,'button','Refresh','save shift_gui');
shifty=add(shifty,'button','Go','auto_shift_DCE_images_callback(params)');
shifty=add(shifty,'button','Quit','close');

shifty=parametergui(shifty);


