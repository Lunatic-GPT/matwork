function show_wash_in_out

shifty=parameter('Show Wash in&out');

shifty=add(shifty,'directoryname','Select image directory',fullfile('.','DCE_temp'));
shifty=add(shifty,'directoryname','Select data directory',fullfile('.','Data'));

shifty=add(shifty,'int','First slice','1');

dir_str = dir('DCE_temp/*_1.mat');
ns = length(dir_str);
if ns == 0
    ns = 3;
end

shifty=add(shifty,'int','Total number of slices',num2str(ns));
shifty=add(shifty,'button','Show','show_wash_in_out_gui(params)');
shifty=add(shifty,'button','Quit','close');

shifty=parametergui(shifty);
