function fm=Daubechies4(m,l)

if ~exist('l','var')
    l=size(m,1);
end
if ndims(m)<3
    
   if l>1
       
        
    for i=1:l
        m(1:l,i)=Daubechies4_1D(m(1:l,i));
    end
    
    for i=1:l
        m(i,1:l)=Daubechies4_1D(m(i,1:l));
    end
   
    fm=Daubechies4(m,l/2);
   else
       fm=m;
    
   end
else
    
    error('not supported yet');
end
    
    

