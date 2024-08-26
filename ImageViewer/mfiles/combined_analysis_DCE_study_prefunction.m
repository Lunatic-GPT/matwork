function combined_analysis_DCE_study_prefunction


analysis=parameter('Detailed DCE analysis');

analysis=add(analysis,'directoryname','Select data directory',fullfile('.','Data'));

analysis=add(analysis,'directoryname','Select image directory',fullfile('.','DCE_temp')); 
analysis=add(analysis,'int','First slice','1');
analysis=add(analysis,'int','Total number of slices','3');

analysis=add(analysis,'bool','use 3d ROI',true);


analysis=add(analysis,'float','right cut-off','36');

analysis=add(analysis,'float','left cut-off','-28.7');
analysis=add(analysis,'button','Analyze!','DCE_cluster_analysis_for_BLAT(params)');
analysis=add(analysis,'button','Quit','close');


parametergui(analysis);
