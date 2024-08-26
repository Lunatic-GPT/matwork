function [V,p1,p2,an,resid]=Path_Orientation_AllPathVox(Path_mFile,dir4allvox)
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

V = zeros([sz,3],'single');

    
for i=1:length(i_ind_path)
 
       iind=i_ind_path{i}(:);
       
   %    p1{i}(j,:)=ind2subb(sz,ind{i}(iind(1)));
   %    p2{i}(j,:)=ind2subb(sz,ind{i}(iind(end)));
       
       sub=ind2subb(sz,ind{i}(iind));
       
       if dir4allvox
          
            [pca_res,score,latent]=pca(sub(:,:));
       end
       
       for k=1:length(iind)
           
           if k<3 || k>length(iind)-2
               continue;
           end
           
           if ~dir4allvox
             [pca_res,score,latent]=pca(sub(k-2:k+2,:));
           end
           
           V(sub(k,1),sub(k,2),sub(k,3),:)= pca_res(:,1);
           
           %%
           p1{i}(k,:)=sub(k,:)+pca_res(:,1)';
           p2{i}(k,:)=sub(k,:)-pca_res(:,1)';
           
           [an{i}(k),phi]=unitVec2thetaPhi(pca_res(:,1).*voxsize(1:3)');
           
           resid{i}(k)=1-latent(1)/sum(latent);
           
           %%
       end

end





