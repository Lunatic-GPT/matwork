function Analyze

analysis=parameter('Analyze');

analysis=add(analysis,'directoryname','Select image directory',fullfile('./','DCE_temp'));
analysis=add(analysis,'directoryname','Select data directory',fullfile('./','Data'));

analysis=add(analysis,'int','First slice','1');

dir_str = dir('DCE_temp/*_1.mat');
ns = length(dir_str);
if ns == 0
    ns = 3;
end
    
analysis=add(analysis,'int','Total number of slices',num2str(ns));
analysis=add(analysis,'bool','Perform motion correction',false);
%analysis=add(analysis,'bool','use 3d ROI',true);

analysis=add(analysis,'float','right cut-off','36');

analysis=add(analysis,'float','left cut-off','0');

analysis=add(analysis,'button','Go','Analyze_callback(params)');

analysis=add(analysis,'button','Quit','close');

analysis=parametergui(analysis);


