function res=roiRing(m,n)
 
% generates a ring with n voxels around clusters in m;



m0=m;
res=0*m0;

for i=1:size(m0,3)
    m=m0(:,:,i);
    for rep=1:n
        
        m=conv2(double(m>0),ones(3,3));
        m=m(2:end-1,2:end-1);
        
    end
    
    res(:,:,i)=m>0&m0(:,:,i)==0;
    
end



    
