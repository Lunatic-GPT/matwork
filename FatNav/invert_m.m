
function res=invert_m(m)

res=0*m;
for i=1:size(m,3)
    
    tmp=m(:,:,i);
   tmp(end+1,:)=[0,0,0,1];
  itmp=inv(tmp);
  res(1:3,:,i)=itmp(1:3,:);

end
