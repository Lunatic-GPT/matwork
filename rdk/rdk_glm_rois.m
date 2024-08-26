function rdk_glm_rois
sid = 'K004';
add_rdk_subject(sid);
if ~strcmp(sid,'K004')
!3dcalc -a 'rdk_glm_REML+orig[6]' -b 'rdk_glm_REML+orig[10]' -c 'rdk_glm_REML+orig[14]' -expr '(a+b+c)/3*226' -prefix glm_REML_coef_coh234
!3dcalc -prefix mask_REML_coh234_thr3.29 -a mask_REML_coh2thr3.294+orig -b mask_REML_coh3thr3.294+orig -c mask_REML_coh4thr3.294+orig -expr 'and(a,b,c)'
!3dTcat glm_REML_coef_coh234+orig mask_REML_coh234_thr3.29+orig. -prefix glm_REML_coef_coh234b
maskcalc('mask_REML_coh234_thr3.29+orig','step(a)',7,'mask_REML_coh234_thr3.29_clst');
else
   %!3dcalc -a 'rdk_glm_REML+orig[1]' -b 'rdk_glm_REML+orig[3]' -c 'rdk_glm_REML+orig[5]' -expr '(a+b+c)/3' -prefix glm_REML_coef_coh234
   !3dcalc -prefix mask_REML_coh234_thr3.29 -a mask_REML_coh1thr3.294+orig -b mask_REML_coh2thr3.294+orig -c mask_REML_coh3thr3.294+orig -expr 'and(a,b,c)'
   !3dTcat glm_REML_coef_coh234+orig mask_REML_coh234_thr3.29+orig. -prefix glm_REML_coef_coh234b
maskcalc('mask_REML_coh234_thr3.29+orig','step(a)',7,'mask_REML_coh234_thr3.29_clst');
end