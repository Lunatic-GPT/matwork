function res=line2roi(p1,p2,sz,rad)

res=zeros(sz);

if length(sz)==2
    p1(3)=1;
    p2(3)=1;
    sz(3)=1;
end

dist=sos(p1(:)-p2(:),1);

x=linspace(p1(1),p2(1),2*ceil(dist));

y=linspace(p1(2),p2(2),2*ceil(dist));
z=linspace(p1(3),p2(3),2*ceil(dist));

rlist=linspace(0,rad,round(rad/0.5)+1);

perp=perp2Line3D(p1-p2);

for ix=1:length(x)
    
    for r=rlist
        philist=linspace(0,2*pi,round(4*pi*r)+2);
        philist=philist(1:end-1);
        for phi=philist
            
            pos=[x(ix),y(ix),z(ix)]+perp(2,:)*sin(phi)*r+perp(1,:)*cos(phi)*r;
            
            i1=round(pos(1));
            i2=round(pos(2));
            i3=round(pos(3));
            if i1<=0||i1>sz(1) ||i2<=0||i2>sz(2)||i3<=0||i3>sz(3)
                continue;
            end
              res(i1,i2,i3)=1;
        end
        
      
    end
end

