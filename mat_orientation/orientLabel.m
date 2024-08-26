
function ori=orientLabel(v)
% v is the rotation matrix;

label='LRPASI';  % the negative side


for i=1:size(v,2)
   
    [tmp,ind]=max(abs(v(:,i)));

    j=ind*2+(sign(v(ind,i))-1)/2;
    ori(i)=label(j);
%     f2(i)=ind;
    
end
