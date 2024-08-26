function dcm2pvs_manualACPC(dpattern,server,step_start)
if ~exist('server','var') || isempty(server)
 server='andrew';
end

if ~exist('step_start','var')
    step_start=1;
end

th=tic;
  dall=name4pat(dpattern);
 dall=str2cell(dall);
 dname=dall{1};
 
 if step_start==1
  if ~exist([dname,'.nii.gz'],'file')
   dcm2nii(dname);
  end
 end
 
     if step_start<=2
         get_ac_pc(dname,1);
     end
     [roi,ind_mid]=get_ac_pc(dname,2);
     
 
 
if step_start<=3
  pvsseg_remote(dname,server);
end

if step_start==4
  pvsseg_remote(dname,server,1); %download only
end

if step_start<=5
    prefix=gen_pvsmask(sprintf('prob_%s.nii.gz',filename(dname)));
    process_pvs([prefix,'.nii.gz']);
else
    prefix=gen_pvsmask(sprintf('prob_%s.nii.gz',filename(dname)));
end

load(fullfile(prefix,[prefix,'_stat.mat']),'pvs');

[icenter,perp]=get_points2surface(roi);

orient=get_orient_from_nii([dname,'.nii.gz']);
n=zeros(1,3);
nperp=zeros(1,3);
norm=zeros(3,3);
center=zeros(3,3);
ipvs_perp=cell(1,3);
sgn=[-1,1,0];
    for lr=1:3
        center(:,lr)=ijk2xyz([ind_mid+20/orient.voxsize(1)*sgn(lr),icenter+perp*15/orient.voxsize(2)],orient);
       
        res(lr)=Slice_P2PVS(pvs,orient,center(:,lr),[0,perp]);
        
        if ~isempty(res(lr).th)
            [nperp(lr),n(lr),norm(:,lr),ipvs_perp{lr}]=find_orient_max_PVS(res(lr),10); %find the slice tilt with the most perpendicular PVS
        end
    end

save(sprintf('PCSlice_%s.mat',dname),'nperp','n','norm','center','ipvs_perp','res');
fprintf('Perpendicular/Total PVSs: %d/%d (r); %d/%d (l); %d/%d (both)\n',nperp(1),n(1),nperp(2),n(2),nperp(3),n(3));

for i=1:3
 save_slicePosition(i,center(:,i),norm(:,i));
end

label='lrb';
for lr=1:3
 outname=sprintf('%s_reslice_%s.mat',dname,label(lr));
 relice_T2_mask(pvs,sprintf('%s.nii.gz',dname),center(:,lr),norm(:,lr),ipvs_perp{lr},outname);
end
disp('');

toc(th);


function save_slicePosition(wip,pos,norm)

if exist('Y:\','dir')
 fname=sprintf('Y:\SlicePosition%d.txt',wip);
else
 fname=sprintf('SlicePosition%d.txt',wip);
end

fid=fopen(fname,'w');

fprintf(fid,'position=(%f,%f,%f)\n',pos);
fprintf(fid,'normal=(%f,%f,%f)\n',norm);
fclose(fid);

copyfile(fname,filename_append(fname,'_base'));


function [nperp_max,n_max,norm,ipvs_perp]=find_orient_max_PVS(res,thr)
% return the index array of intersecting PVSs and their tilt angles wrt the imaging plane

for i=1:length(res.angle)
    nperp(i)=sum(res.angle{i}<thr);
    n(i)=length(res.angle{i});
end
[nperp_max,imax]=max(nperp);
norm = res.norm(:,imax);
ipvs_perp=res.ipvs{imax}(res.angle{imax}<thr);
Norm2SliceOrient(norm);
n_max=n(imax);

    
function [center,norm]=get_points2surface(roi)

a=clusterize2(roi);

pos_in=get_roi_endpoints(a==1);

tang=diff(pos_in,1);
norm=[-tang(2),tang(1)]/sqrt(tang(2)*tang(2)+tang(1)*tang(1));
center = round(mean(pos_in,1));


function pos=get_roi_endpoints(roi)

sub=ind2subb(size(roi),find(roi>0));

[~,xmin]=min(sub(:,1));
[~,xmax]=max(sub(:,1));

pos=[sub(xmin,:);sub(xmax,:)];


%%




