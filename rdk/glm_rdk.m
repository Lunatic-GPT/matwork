p = parameter('GLM analysis in afni');

str = 's1_dtr_norm+orig';
for i =2:10
  str = sprintf('%s,s%d_dtr_norm+orig',str,i);
end 

p = add(p,'string','time series names',fnames);

p = add(p,'filename','mask','brainmask+orig');
p = add(p,'filename','stim_file',[]);

p = add(p,'string','general linear tests','');
p = add(p,'button','browse','browse_callback(params,''general linear tests'')');

p = add(p,'bool','add global signal?',false);
p = add(p,'float','1D TR',2);
p = add(p,'int','baseline order',3);
p = add(p,'string','output prefix','glm');

p = add(p,'button','Save','disp(''parameters saved to params_glm_afni_gui.mat'');save(''params_glm_afni_gui.mat'',''params'');');
p = add(p,'button','Load','load_callback(params)');

p = add(p,'button','Calculate','glm_afni(params);');
p = add(p,'button','1D calculate','glm_afni1D(params);');
p = add(p,'button','Close','close');