function Analyze_callback(params)

set(gcbo,'Enable','off');
auto_shift_DCE_images_callback(params);

DCE_cluster_analysis_for_BLAT(params);

set(gcbo,'Enable','on');