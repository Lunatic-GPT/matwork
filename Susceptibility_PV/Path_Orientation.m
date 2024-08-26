function [p1,p2,an,resid]=Path_Orientation(Path_mFile,lseg)
% p1: 1*npath cell; each cell element is a n*3 matrix for n [x,y,z]
% coordinates
% p2: same as p1, the other end of the path segment
% an: angle of the segment as determined from the fit
% resid: the residual for fitting a straight line to the path segment
% path_dir: folder containing exported pathes
% len: divide the paths into segments of length len


i_ind_path=ri(Path_mFile,'','','i_ind_path');
voxsize=ri(Path_mFile,'','','voxsize');

ind=ri(Path_mFile,'','','ind');
c=ri(Path_mFile,'','','c');
sz=size(c);

m=c*0;

    
for i=1:length(i_ind_path)

   [l,lall]=pathLength(i_ind_path{i},ind{i},voxsize,sz);

   if ~isinf(lseg)
       nseg=floor(l/lseg);
   else
       nseg=1;
   end
   
   i_i_ind_path = pair_vox2path(ind{i},sz,i_ind_path{i});
    
   for j=1:nseg
       
       subind=find(lall>=lseg*(j-1)&lall<=lseg*j);
       
       iind=i_ind_path{i}(subind);
       
   %    p1{i}(j,:)=ind2subb(sz,ind{i}(iind(1)));
   %    p2{i}(j,:)=ind2subb(sz,ind{i}(iind(end)));
       
       for k=1:length(subind)
           iind=[iind,find(subind(k)==i_i_ind_path)];
       end
       iind=unique(iind);
       
       m(ind{i}(iind))=j;
       sub=ind2subb(sz,ind{i}(iind));
       
       
       
      [pca_res,score,latent]=pca(sub);
      
      
      % find the voxel in cluster that is nearest to mean(sub,1)

      sub_c=ind2subb(sz,ind{i});
      
      dist=sos(repmat(mean(sub,1),[size(sub_c,1),1])-sub_c,2);
      
      [tmp,imin]=min(dist);
       
      p1{i}(j,:)=sub_c(imin,:)+pca_res(:,1)';
      p2{i}(j,:)=sub_c(imin,:)-pca_res(:,1)';
      
      [an{i}(j),phi]=unitVec2thetaPhi(pca_res(:,1).*voxsize(1:3)');
       
      resid{i}(j)=1-latent(1)/sum(latent);
       
   end
   
end

prefix=strtok(Path_mFile,'.');

save([prefix,'_seg.mat'],'m','voxsize');





