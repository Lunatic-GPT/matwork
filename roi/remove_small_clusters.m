function remove_small_clusters(fpat,thr,twoD)

%fpat: file pattern
%thr: minimum number of voxels in the cluster to survive.
%twoD: form clusters slice by slice or in 3D.

a=dir(fpat);

for i=1:length(a)
    
roi=load_untouch_niigz(a(i).name);

if twoD
wm=clusterize2_2d(roi.img>0,thr);
else
    wm=clusterize2_2d(roi.img>0,thr);
end
roi.img=wm;
prefix=strtok(a(i).name,'.');
save_untouch_niigz(roi,[prefix,'_clustered']);

end

 
 