function shift_DCE_images_stand_alone_prefunction
% SHIFT_DCE_IMAGES_STAND_ALONE_PREFUNCTION lets you denote the shifts in
% the DCE series and then actually perform the shift

%load root_directory
shifty=parameter('DCE image shifting (x/y shifts)');

shifty=add(shifty,'directoryname','Select image directory name',pwd);
%shifty=setcallback(shifty,'Select image directory name', 'save shift_gui');


shifty=add(shifty,'int','2-1x','0');
%shifty=setcallback(shifty,'2-1x','save shift_gui');
shifty=add(shifty,'int','2-1y','0');
%shifty=setcallback(shifty,'2-1y','save shift_gui');
shifty=add(shifty,'int','2-3x','0');
%shifty=setcallback(shifty,'2-3x','save shift_gui');
shifty=add(shifty,'int','2-3y','0');
%shifty=setcallback(shifty,'2-3y','save shift_gui');
shifty=add(shifty,'int','2-4x','0');
%shifty=setcallback(shifty,'2-4x','save shift_gui');
shifty=add(shifty,'int','2-4y','0');
%shifty=setcallback(shifty,'2-4y','save shift_gui');
shifty=add(shifty,'int','2-5x','0');
%shifty=setcallback(shifty,'2-5x','save shift_gui');
shifty=add(shifty,'int','2-5y','0');
%shifty=setcallback(shifty,'2-5y','save shift_gui');
shifty=add(shifty,'int','2-6x','0');
%shifty=setcallback(shifty,'2-6x','save shift_gui');
shifty=add(shifty,'int','2-6y','0');
%shifty=setcallback(shifty,'2-6y','save shift_gui');

shifty=add(shifty,'string','Enter the subject ID','N84432');
%shifty=setcallback(shifty,'Enter the subject ID','save shift_gui');

shifty=add(shifty,'int','Total number of slices','3');
%shifty=setcallback(shifty,'Total number of slices', 'save shift_gui');

shifty=add(shifty,'button','Refresh','save shift_gui');
shifty=add(shifty,'button','Shift!','shift_DCE_images_stand_alone(params)');
shifty=parametergui(shifty);


