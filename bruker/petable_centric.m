function  t2=petable_centric(N,nr,tbfactor)
% t2=petable_centric(N,nr,tbfactor)
if ~exist('nr','var')
    nr=1;
end

if ~exist('tbfactor','var')
    tbfactor=1;
end


petable=-N/2:N/2-1;

[tmp,ind]=sort(abs(petable),'descend');
    
t2=petable(ind);

for i=1:N/tbfactor
    
   jj=(i-1)*tbfactor+1:i*tbfactor;
   
   ttmp=t2(jj);
   
   [tmp,ind]=sort(abs(ttmp),'ascend');
   
   t2(jj)=ttmp(ind);
    
end


t2=repmat(t2,[1,nr]);

if nargout==0
 save_mat_int(t2',sprintf('petable_centric_%d_NR%d',N,nr));
end


