function res=roiCOMNeighbor(m,n)
 
% generates a ring with n voxels around clusters in m;



res=0*m;

l1=floor(n/2);
l2=floor((n-1)/2);


for i=1:size(m,3)
ind=ind2subb(size(m(:,:,1)),find(m(:,:,i)>0));

if ~isempty(ind)
com=round(mean(ind,1));  % center of mass


res(com(1)-l1:com(1)+l2,com(2)-l1:com(2)+l2,i)=1;
end

end

res(m>0)=0;