function [coef,err]=roi_glm(tseries,xmat,mask)
% [coef,err]=roi_glm(tseries,xmat,mask)
% do not include the baseline model in the design matrix.
ts = [];
tseries = str2cell(tseries);
rname = '';
fst = 0;
for i=1:length(tseries)
    [n,tmp] = nv_roi(tseries{i},mask);
    rname = sprintf('%s %d',rname,fst);
    fst = fst + length(tmp);
    ts = [ts;tmp(:)];
end

save ts_roi_glm.1D ts -ASCII 
p = parameter('GLM analysis in afni of 1D file');

p = add(p,'filename','time series names','ts_roi_glm.1D');

p = add(p,'filename','stim_file',xmat);

p = add(p,'string','-concat ',sprintf('''1D: %s''',rname));
p = add(p,'string','general linear tests','');

p = add(p,'bool','add global signal?',false);
p = add(p,'float','1D TR',2);
p = add(p,'int','baseline order',3);
p = add(p,'string','output prefix','results_roi_glm');
glm_afni1D(p);

d= textread('results_roi_glm_REML.1D','','commentstyle','shell');

coef = d(2:2:end-1);
f = d(3:2:end);
err = coef./sqrt(f);





