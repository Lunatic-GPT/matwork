function xmat_itk2dcm(fname)
a=load(fname);
a=reshape(a,[3,4]);
a=[a;0,0,0,1];
m=eye(4);
m(1,1)=-1;
m(2,2)=-1;
a=m*a*m;

save(filename_append(fname,'_new'),'a','-ascii');



