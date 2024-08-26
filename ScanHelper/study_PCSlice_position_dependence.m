function study_PCSlice_position_dependence
addpath(fullfile(topvsana,'recon_tse'));
sid=get_exclude_sub_motion(1);

i=5;
if sid(i)==0
    disp('excluded');
    return;
end
f_pvs=sprintf('../pvsseg/mask_koji_ver2/pvsMask_PVS%02d.nii.gz',i);
f_T2=sprintf('../T2nii/T2_PVS%02d',i);

[roi,ind_mid]=get_ac_pc(f_T2);


[icenter,perp]=get_points2surface(roi);
pvs=process_pvs(f_pvs,1,0,0);
orient=get_orient_from_nii([f_T2,'.nii.gz']);
n=zeros(1,2);
nperp=zeros(1,2);
norm=zeros(3,2);
center=zeros(3,2);
ipvs_perp=cell(1,2);
sgn=[-1,1];

    for lr=1:2
        center(:,lr)=ijk2xyz([ind_mid+20/orient.voxsize(1)*sgn(lr),icenter+perp*15/orient.voxsize(2)],orient);
        res(lr)=Slice_P2PVS(pvs,orient,center(:,lr),[0,perp]);
        
        if ~isempty(res(lr).th)
            [nperp(lr),n(lr),norm(:,lr),ipvs_perp{lr}]=find_orient_max_PVS(res(lr),10); %find the slice tilt with the most perpendicular PVS
        end
    end


save(sprintf('study_PCSlice_PVS%02d',i),'nperp','n','norm','res','ipvs_perp');
label='rl';
for lr=1:2
[~,ind]=max(nperp(:,lr));
outname=sprintf('%s_%s.mat',filename(f_T2),label(lr));
relice_T2_mask(pvs,sprintf('%s.nii.gz',f_T2),center(:,lr),norm(:,lr),ipvs_perp{lr},outname);
end
disp('');

function [nperp_max,n_max,norm,ipvs_perp]=find_orient_max_PVS(res,thr)
% return the index array of intersecting PVSs and their tilt angles wrt the imaging plane

for i=1:length(res.angle)
    nperp(i)=sum(res.angle{i}<thr);
    n(i)=length(res.angle{i});
end
[nperp_max,imax]=max(nperp);
norm = res.norm(:,imax);
ipvs_perp=res.ipvs{imax}(res.angle{imax}<thr);

n_max=n(imax);
    
function [center,norm]=get_points2surface(roi)

a=clusterize2(roi);

pos1=get_roi_endpoints(a==1);
pos2=get_roi_endpoints(a==2);

pos_in=pos1;
pos_out=pos2;

if mean(pos1(:,2),1)>mean(pos2(:,2),1)
   pos_in=pos2;    
   pos_out=pos1;
end

tang=diff(pos_in,1);

norm=[-tang(2),tang(1)]/sqrt(tang(2)*tang(2)+tang(1)*tang(1));

 
center = round(mean(pos_in,1));

%xbound=ceil((size(roi,2)-mean(pos_in(:,2),1))/perp(2));
%{
for i=1:xbound
   res(i,:) = round(mean(pos_in,1)+i*perp);
  v1=res(i,:)-pos_out(1,:);
  v2=pos_out(2,:)-pos_out(1,:);
  crs=cross([v1,0],[v2,0]);
  if crs(3)<0
      res(end,:)=[];
    break;
    
  end
 
end

res=unique(res,'rows');
%}
function pos=get_roi_endpoints(roi)

sub=ind2subb(size(roi),find(roi>0));

[~,xmin]=min(sub(:,1));
[~,xmax]=max(sub(:,1));

pos=[sub(xmin,:);sub(xmax,:)];


%%




