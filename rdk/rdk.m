function rdk
p = parameter('rdk analysis');
p = add(p,'button','go to subject','edit add_rdk_subject');

p = add(p,'button','reorder dicom','edit reorder_dicom_gui');

p = add(p,'button','sort time series','edit ts_sort_bin_gui');
p = add(p,'button','correlation analysis','edit rdk_cc_gui');
p = add(p,'button','Retinotopy/hMT localizer','edit rdk_retino_loc');
p = add(p,'button','glm analysis','glm_afni_gui');
p = add(p,'button','Inspect sorted time course in roi','edit plot_sorted_ts');
p = add(p,'button','Calcuate roi mean time course','edit batch_plot_single'); 
p = add(p,'button','File Rename','edit rdk_glm_rois'); 
p = add(p,'button','Plot results','edit plot_group'); 
p = add(p,'button','Quit','close');
parametergui(p);