function m=neighbors_mask(pos,dist,sz)
%m=neighbors_mask(pos,dist,sz)
%pos are one-based integers; n*3
% dist 

th=(0:19)*2*pi/20;

dpos=diff(pos,1,1);
m=zeros(sz);
for i=1:size(pos,1)-1
   n1=dpos(i,:);
 [tmp,t]=max(abs(n1));
 
 j=1:3;
 j(t)=[];
 n2=ones(1,3);
  
 n2(t)=-sum(n1(j).*[1,1])/n1(t);
 
 n2=n2/sqrt(sum(n2.^2)); 
 n1=n1/sqrt(sum(n1.^2));
 n3=cross(n1,n2);
    
 orig=pos(i:i+1,:);
 for k=1:length(th)
    
     for l=1:size(orig,1)
      coord=round(orig(l,:)+dist*(n2*sin(th(k))+n3*cos(th(k))));
      
      coord(coord<1)=1;
      coord(coord>sz)=sz(coord>sz);
      
      m(coord(1),coord(2),coord(3))=1;
      
     end
 end
 
 
end