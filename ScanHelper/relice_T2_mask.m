function relice_T2_mask(pvs,f_T2,center,norm,ipvs_perp,outname)
% center: in mm; 
% norm: in dicom convention

T2=do_SliceT2_T2C_parietal(center,norm,f_T2);  %debug
orient=get_orient_from_nii(f_T2);
o_ref=get_o_ref(norm,center);
m=find_intercept2(pvs,o_ref,orient);

m_perp=remove_tilt_pvs(m,ipvs_perp);

save(outname, 'T2', 'm', 'm_perp');

function c2=find_intercept2(pvs,o_ref,o_ind)
% norm: 1*3
% norm_seg: 3*n cells
% the coordinate system assumes rotmat=eye(3) and pos =[0,0,0] for the T2
% images;


sz_interp=[512,512,1];
c2=zeros(sz_interp);

ipvs=[];
for i=1:length(pvs)
    
    xyz=ijk2xyz(pvs(i).subind(1:pvs(i).nvox_path,:),o_ind);
    res=dist2plane(xyz,o_ref.center,o_ref.rotmat(:,3));
    
    iint=find(sign(res(1:end-1)).*sign(res(2:end))<=0);
    if length(iint)>1 || isempty(iint)
        continue;
    end
    
    if abs(res(iint))<abs(res(iint+1))
        %an(count)=angle_bw_2vec(norm_seg{i}(:,iint)',norm);      
        iint2=iint;
    else
        %an(count)=angle_bw_2vec(norm_seg{i}(:,iint+1)',norm);      
        iint2=iint+1;
    end
    ijk=xyz2ijk(xyz(iint2,:),o_ref);
        
    
    c2(ijk(1),ijk(2))=i;
 %   T2b(ijk(1),c2(ijk(2)))=T2(sub(iint2,1),sub(iint2,2),sub(iint2,3));
        ipvs(end+1)=i;
end

%disp(ipvs);

function m=remove_tilt_pvs(m,ipvs_perp)

ipvs=unique(m(:));
ipvs(ipvs==0)=[];
ipvs(isnan(ipvs))=[];

ipvs_tilt=setdiff(ipvs,ipvs_perp);
if ~isempty(ipvs_tilt)
    for i=ipvs_tilt'
        m(m==i)=0;
    end
end


function o_ref=get_o_ref(norm,center)
     rotmat=norm2rotmat(norm);
    
    o_ref.voxsize=[0.4,0.4,1];
    o_ref.center=center;
    
    o_ref.pos=center2pos(o_ref.voxsize,rotmat,[512,512,1],center);
    o_ref.rotmat=rotmat;
    o_ref.oreint='LPI';
    
function T2b=do_SliceT2_T2C_parietal(center,norm,f_T2) % for debug 

    T2=ri(f_T2);
    o_T2=get_orient_from_nii(f_T2);
     rotmat=norm2rotmat(norm);
    o_ref.voxsize=[0.4,0.4,1];
    o_ref.center=center;    
    o_ref.pos=center2pos(o_ref.voxsize,rotmat,[512,512,1],center);
    o_ref.rotmat=rotmat;
    o_ref.oreint='LPI';  
    T2b=reslice_with_orient([512,512,1],o_ref,o_T2,T2);

    
function res=norm2rotmat(norm)

c=norm(:);
a=[1,0,0]';
a=a-c*sum(c.*a);
a=a/sos(a,1);
b=cross(c,a,1);
res=[a,b,c];

function res=dist2plane(v,center,norm)
% v=n*3;
% center: 1*3;
% norm=1*3

res=sum((v-center(:)').*norm(:)',2);

   