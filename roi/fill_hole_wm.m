function fill_hole_wm(fpat,thr)

a=dir(fpat);

for i=1:length(a)
    
roi=load_untouch_niigz(a(i).name);
wm=clusterize2_hole(roi.img==2,thr);

roi.img=wm;
prefix=strtok(a(i).name,'.');
save_untouch_niigz(roi,[prefix,'_nohole']);

end

 
 