function prob2mask_nii(fpattern,thr)

if ~exist('thr','var')
    thr=0.8;
end

d=name4pat(fpattern);
for i=1:length(d)
    fname=filename(d{i});
    prefix=strtok(fname,'.');
    fprintf('%s\n',prefix);  
    [id,suf]=strtok(prefix,'_');
    if exist(['mask_',suf(2:end),'.nii.gz'],'file')
        continue;
    end
    
    
    a=load_untouch_niigz(d{i});
    
    a.img=int16(a.img>thr);
    n=clusterize2(a.img>0);
    fprintf('number of clusters is %d\n',max(n(:)));
    save_untouch_niigz(a,['mask_',suf(2:end)]);
    
end