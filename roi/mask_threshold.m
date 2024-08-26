function roi=mask_threshold(d,m,sd_scale,csize,nrep,pos_neg)
%  roi=mask_threshold(d,m,sd_scale,csize,nrep,pos_neg)
% generate mask based on thresholding
% d: the data
% m: only voxels inside m will be considered
% sd_scale: thr=sd_scale*std;  CI: 95%: sd_scale = 1.96 (two tail); 1.64 (one tail)
% nrep: std calculated nrep times excluding voxels outside the treshold
% (non-background)
% pos_neg: -1, mask hypo; 1 mask hyper; 0 both 

x=d(m>0);

sig=Inf;
mn=0;
disp('Thresholds:');
for i=1:nrep
    
    nsd=sig*sd_scale;
    m2=m&abs(d-mn)<=nsd;
    m3=clusterize2(m2,csize);
    mn=mean(d(m3>0));
    sig=std(d(m3>0));
    disp(sig*sd_scale);
end


if pos_neg==-1
    roi=m>0&d<=mn-sig*sd_scale;
elseif pos_neg==1
    roi=m>0&d>=mn+sig*sd_scale;
else
    roi=m>0&abs(d-mn)>=sig*sd_scale;
end
roi=clusterize2(roi,csize);
