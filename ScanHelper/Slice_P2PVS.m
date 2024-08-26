function res=Slice_P2PVS(pvs,orient,center,norm_ave)
%  res=Slice_P2PVS(pvs,orient,center)
% pvs: structure of PVSs
% orient: orient struct of the PVS scan
% center: center of the PC slice
% norm_ave: the normal direction around which a small 3D angular region is defined to search PVS 
% 10/15/2020: add th=0 and phi =0 to the search.


d2r=pi/180;

th=(0:10)*d2r;
th_all=[];
ph_all=[];

count=0;
for i=1:length(th)
    dth=th(2)-th(1);
   
    dphi=dth/sin(th(i));
    
    phi=0:dphi:2*pi;
    
    for j=1:length(phi)    
       norm=get_norm(norm_ave,th(i)/d2r,phi(j)/d2r);   
       count=count+1;
       [ag{count},ipvs{count}]=find_PVS_intersect_slice(pvs,orient,center,norm);               
     
       norm_all(:,count)=norm; 
       th_all(count)=th(i)/d2r;
       ph_all(count)=phi(j)/d2r;
    end
      
end

res.th=th_all;
res.ph=ph_all;
res.ipvs=ipvs;
res.angle=ag;
res.norm_ave=norm_ave;
res.norm=norm_all;


function res = get_norm(norm_ave,th,phi)
        [th0,phi0]=unitVec2thetaPhi(norm_ave(:));
        m=transform_matrix_rotation(th0,phi0);
        
        vec=thetaPhi2unitVec(th,phi);        
        res=m*vec';
        








