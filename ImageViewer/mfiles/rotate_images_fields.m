function d2=rotate_images_fields(d,ax,phi)
% d2=rotate_images(d,ax,phi)
% ax: 1,2, or 3
% ph: angle in degrees

an=phi*pi/180;

c=cos(an);
s=sin(an);

ind=1:3;

ind(ax)=[];

mat=eye(3);
mat(ind(1),ind(2))=s;
mat(ind(2),ind(1))=-s;
mat(ind(1),ind(1))=c;
mat(ind(2),ind(2))=c;


d2=d*0;
sz=size(d);


for i=1:size(d,1)
    tic;
    for j=1:size(d,2)
     %   disp([i,j]);
        for k=1:size(d,3)
            
            pos=ijk2pos([i,j,k],sz(1:3));  % coordinate after rotation
            
            xyz_new=pos2ijk(mat*pos',sz(1:3));  % xp is the coordinates in the original data space
                  
             xyz_l=floor(xyz_new);
             xyz_u=xyz_l+1;
             
             xyz=[xyz_l(:),xyz_u(:)];
             if any(xyz_l<1)||any(xyz_u>sz(1:3))
                 continue;
             end
             
             for l=1:3
              d2(i,j,k,l)=trilinear_interpolation(xyz_new,xyz,d(xyz_l(1):xyz_u(1),xyz_l(2):xyz_u(2),xyz_l(3):xyz_u(3),l));
            
             end
             
             d2(i,j,k,:)=mat\vec(d2(i,j,k,:));
             
        end
    end
    disp([i,toc]);
end

           
function ijk=pos2ijk(pos,N)

  dv=2./N;
  ijk=(pos'+1-dv/2)./dv;
  

function pos=ijk2pos(ijk,N)

dv=2./N;
pos=-1+dv/2+dv.*ijk;






