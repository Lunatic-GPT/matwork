function pvs=PVS_quantify(pvs,nii,prefix)
% 
npath=length(pvs);
v_allseg=[];  % v for all segments
  
voxsize=nii.hdr.dime.pixdim(2:4);
orient=get_orient_from_nii(nii);
sz=size(nii.img);

d_allseg=[]; %diameter for all segments.
lthr=0.8; %
diam=zeros(sz);
dir_seg=zeros([sz,3]);
for i=1:npath
 
    posc=pvs(i).subind;  
    posc=posc.*repmat(voxsize(:)',[size(posc,1),1]);
     
    pospath=posc(1:pvs(i).nvox_path,:);

    v=length(pvs(i).ind)*prod(voxsize);
   
    lnorm=zeros(1,pvs(i).nvox_path);
    
    l=0;
 
    for j=1:pvs(i).nvox_path-1
        l=l+sqrt(sum((pospath(j,:)-pospath(j+1,:)).^2));
        lnorm(j+1)=l;   
    end
    pvs(i).l=l; 
    pvs(i).lnorm=lnorm/l;
     pvs(i).v=v;
     
    l=l+sos(pospath(1,:)-pospath(2,:))/2+sos(pospath(end-1,:)-pospath(end,:))/2;
    r=sqrt(v/l/pi);
    pvs(i).r=r;
    
      lseg=zeros(1,pvs(i).nvox_path);
      for j=2:pvs(i).nvox_path-1
         
          lseg(j)=0.5*(sos(pospath(j,:)-pospath(j-1,:)))+0.5*(sos(pospath(j,:)-pospath(j+1,:)));
          
      end
      lseg(1)=sos(pospath(2,:)-pospath(1,:));
      lseg(end)=sos(pospath(end,:)-pospath(end-1,:));
      
      pvs(i).lseg=lseg;
      pvs(i).tortuosity=tortuosity(pospath(1:pvs(i).nvox_path,:),voxsize);
      
      vseg=zeros(1,pvs(i).nvox_path);
          
      pvs(i).ind_pathVoxelLabel = 0*pvs(i).ind;
      
      for j=1:length(pvs(i).ind)
            
          dist=sos(posc(j,:)-pospath); 
          [mdist,ind_min]=min(dist);
          vseg(mdist==dist)=vseg(mdist==dist)+prod(voxsize)/sum(mdist==dist); 
          pvs(i).ind_pathVoxelLabel(j)=ind_min;
        
      end
      
     % fill the dir_seg image 
       for j=1:length(pvs(i).ind)
         si=pvs(i).subind(j,:);
         ipathv=pvs(i).ind_pathVoxelLabel(j);
         dir_seg(si(1),si(2),si(3),:)= pvs(i).dir_seg(:,ipathv);
       end
       
      pvs(i).dpath=sqrt(vseg./lseg/pi)*2;
      
       vsegclust=0*pvs(i).ind;
       for j=1:length(pvs(i).ind)
          dist=sos(posc(j,:)-pospath);
          
          mdist=min(dist);
          vsegclust(j)=mean(vseg(mdist==dist));
          diam(pvs(i).ind(j))=mean(pvs(i).dpath(mdist==dist));
       end
      
      v_allseg(end+1:end+length(vseg))=vseg;
       
      pvs(i).vseg=vseg;  %volume of each segment
      pvs(i).vclust=vsegclust; % volume of the cluster the voxel is associated with
      
      
      if l>=lthr
       d_allseg=[d_allseg,pvs(i).dpath];
      end
      
  
end
% only l, v , and dpath are used in pvs_stat
save([prefix,'_stat.mat'],'pvs','lthr','d_allseg','v_allseg','diam','dir_seg');
if ~isempty(orient)
  save([prefix,'_stat.mat'],'-append','-struct','orient');
end

%% save orientation

img=zeros([sz(1:3),1]);
%image 1: pvs index,
%2:4 - x,y,z compoents of the mean dir vecotr of the entire PVS.
%5:7 - x,y,z compoents of the dir vector of each segment; dicom convention.
dir_whole=zeros([sz(1:3),3]);
dir_seg=zeros([sz(1:3),3]);
for i=1:length(pvs)

    ind=pvs(i).ind;
    norm=pvs(i).norm;
    img(ind)=i;
    dir_whole=setv_ind(dir_whole,ind,norm);

    for j=1:pvs(i).nvox_path
        
       ind4VoP=pvs(i).ind(pvs(i).ind_pathVoxelLabel==j);
       
       dir_seg=setv_ind(dir_seg,ind4VoP,pvs(i).dir_seg(:,j));
     %  time_left(j,pvs(i).nvox_path, toc);
    end           
end
nii.img=cat(4,img,dir_whole,dir_seg);
save_untouch_niigz(nii,[prefix,'_PVSDir']);
            

function [Im]=sos(iX)

% room sum square of the last non-single dim.

if ndims(iX)>2

   Im=sqrt(sum((abs(iX).^2),ndims));

else
    
    if size(iX,2)==1
     Im=sqrt(sum((abs(iX).^2),1));
    else
    Im=sqrt(sum((abs(iX).^2),2));
    end
end

