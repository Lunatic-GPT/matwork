function c=centerFinder(cd,rad,range,circ_square)
% cd: complex data
% rad: radius of circle mask or width of a square (2*rad+1).
% range: the size of the square around center to search the center is
% (2*range+1)*(2*range+1)
% if circ_square = 0: use circle otherwise use a square.

c0=ceil((size(cd)+1)/2);

rl=zeros(2*range+1,2*range+1);

for i=-range:range
    for j=-range:range
             
    if circ_square ==0   
     m=mask_circle(size(cd),rad,c0+[i,j],1);
    else
     m=0*cd;
     m(c0(1)-rad+i:c0(1)+rad+i,c0(2)-rad+j:c0(2)+rad+j)=1;
    end
     rl(i+range+1,j+range+1)=sum(cd(m>0));
     
    end
end


[tmp,ind_min]=min(real(rl(:)));

cshift=ind2subb(size(rl),ind_min);

c=c0+cshift+[-range-1,-range-1];










