function B1map_th2th4th(dlist,theta)
% dlist should contain 3 file names which corresponds to
% theta, 2 theta, and 4 theta respectively or matrix with 3 entries 
% corresponding to theta, 2theta, and 4 theta.
if iscell(dlist)
a= B1map_th2th(dlist(1:2),theta);

b= B1map_th2th(dlist(2:3),2*theta);

thr=90*90/theta/4;

b1=a;
b1(a<thr)=b(a<thr);

save(sprintf('B1map_th2th4th_%s',dlist{1}),'b1');

else
   
    d=dlist;
    if ndims(d)==3
        d=reshape(d,[size(d,1),size(d,2),1,size(d,4)]);
    end
    a= B1map_th2th(d(:,:,:,1:2),theta);

b= B1map_th2th(d(:,:,:,2:3),2*theta);

thr=90*90/theta/4;

b1=a;
b1(a<thr)=b(a<thr);

dlist=strtok(dlist,'.');
save('B1map_th2th4th','b1');

end