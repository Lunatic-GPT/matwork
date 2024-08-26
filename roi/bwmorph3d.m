function m2=bwmorph3d(m,operation,varargin)
m2=0*m;

tmp=sum(sum(m,1),2);

ind=find(tmp(:)>0);

for i=ind'
    
        m2(:,:,i)=bwmorph(m(:,:,i),operation,varargin{:});
        
end


