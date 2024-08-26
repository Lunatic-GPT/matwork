function vroi = sroi2v2(f3dname,sroi_names,subject,lr,file_save,suma,sv,use_norm)
% vroi = sroi2v(f3dname,sroi,subject,lr,file_save[,suma,sv,use_norm])
% This function maps rois defined on the surface to 3d volume.           
% Because of the non-unique nature of the surface to volume mapping. 
% voxels can be classified into the following five types:
%
%     5. voxles that map to a single ROI
%     4. voxels that map to a single ROI and outside any ROI.
%     3. voxels that map to multiple ROIs.
%     2. voxels that map to surface node outside any surface ROIs
%     1. voxels that do not map to any surface node
%
% f3dname: file name of a 3d dataset with one sub-brick.
%          All the voxels with values greater than 0 will be considered
% sroi_name:  1D file names of surface rois (*.lh.1D.roi or *.rh.1D.roi)
% subject: Name of the subject
% lr:      'lh' or 'rh'
% file_save: whether to save the vroi. default: faulse.
% suma: suma directory. Default: '../../../segmentation/<subject>/SUMA'
% sv:  surface volume aligned with the experiment. 
%      Default: '<subject>_SurfVol_ns_Alnd_Exp+orig'
% use_norm: whether to use normals to Surf_A to generate the segments.
%      Default: true;
%
% vroi is an array with the following sub-bricks
% 1. mapping properties of a voxel.  The values are from 1-5 as specified
%    above.   This can be used for thresholding in afni.
% 2. roi index.  non-zero if mapped to only one roi.
% 3. number of mapped rois.
% 4. 1+length(srois) subbricks specifying whether the voxel maps to each of 
%    the rois.  
% 5. 1+length(srois) subbricks specifying the number of 
%    surface nodes within each of the rois mapped to that voxel
% 6. 1+length(srois) subbricks specifying wheather each of the roi has the 
%    maximum number of nodes among the mapped rois for that voxel.
%    This can be used to select the most probable roi among all the rois 
%    that maps to that voxels.

% Note: This function only process one hemisphere. 
%      
%
% History: 1/29/2009, created by X.Z.
% use the afni 3dSurf2Vol mapping directly instead of 3dVol2Surf as in sroiv.  

if ~exist('file_save','var')
    file_save = false;
end

if ~exist('suma','var')
    suma = sprintf('../../segmentation/SUMA');
end

if ~exist('sv','var')
    sv = sprintf('%s_SurfVol_ns_Alnd_Exp+orig',subject);
end

if ~exist('use_norm','var')
   use_norm = true;    
end

       nrois = length(sroi_names);
                       
            [vdata,info] = BrikLoad(f3dname);
            nf = size(vdata,1);
            np = size(vdata,2);
            ns = size(vdata,3);
            vroi = zeros(nf,np,ns,3+3*nrois);
            nn_vdata = zeros(nf,np,ns,nrois+1);   % nrois+1 is used for compatibility with sroi2v   
            f1D = '1Dfile';
            fcount = 'sroi2v2_count';
            
       for i= 1:nrois

            if exist([f1D,'.1D.dset'],'file')
                delete([f1D,'.1D.dset']);
            end
            
            if exist([fcount,'+orig.HEAD'],'file')
                delete([fcount,'+orig.*']);
            end
            
            cmd = sprintf('ROI2dataset -prefix %s -of 1D -input %s',f1D,sroi_names{i});
            unix(cmd);
            if ~use_norm
            cmd = sprintf(['3dSurf2Vol' ... 
             ' -spec %s/%s_both.spec' ...
             ' -surf_A %s.smoothwm.asc'  ...
             ' -surf_B %s.pial.asc'  ...
             ' -sv %s' ...
             ' -grid_parent %s' ...  
             ' -map_func count' ...
             ' -f_steps 10' ...
             ' -sdata_1D  %s.1D.dset' ...
             ' -prefix %s'], ... 
             suma,subject,lr,lr,sv,f3dname,f1D,fcount);
            else
                cmd = sprintf(['3dSurf2Vol' ... 
             ' -spec %s/%s_both.spec' ...
             ' -surf_A %s.smoothwm.asc'  ...
             ' -sv %s' ...
             ' -f_p1_mm 1.0' ...
             ' -grid_parent %s' ...  
             ' -map_func count' ...
             ' -f_steps 10' ...
             ' -sdata_1D  %s.1D.dset' ...
             ' -prefix %s'], ... 
             suma,subject,lr,sv,f3dname,f1D,fcount);
            end
            
             disp('Converting from surface to volume');   
             unix(cmd);
             nn_vdata(:,:,:,i) = BrikLoad([fcount,'+orig']);
        end 
            
      for i1=1:nf
          for i2=1:np
            for i3=1:ns
                             
                nnroi= zeros(1,nrois+1); % number of nodes in each roi that maps to that voxel.  
                for j=1:nrois
                    nnroi(j) = nn_vdata(i1,i2,i3,j);
                end   
                
                nr = 0;
                nmax = max(nnroi);
                for iri = 1:nrois+1
                   if nnroi(iri)>0
                       vroi(i1,i2,i3,3+iri) =1;
                       vroi(i1,i2,i3,3+iri+nrois+1) = nnroi(iri);
                       if nmax == nnroi(iri)
                         vroi(i1,i2,i3,3+iri+2*(nrois+1)) = 1;
                       end
                       nr = nr+1;
                   end
                 end
                 
                vroi(i1,i2,i3,3) = nr; % number of rois.
                                          
                if nr == 1 && nnroi(1+nrois) == 0  %map to a single roi
                    vroi(i1,i2,i3,1) = 5;
                    vroi(i1,i2,i3,2) = find(nnroi(:)>0);
                elseif nr ==2 && nnroi(1+nrois)>=1  % this case never happen here
                    vroi(i1,i2,i3,1) = 4;
                elseif nr > 1  % map to more than one rois
                    vroi(i1,i2,i3,1) = 3;
                elseif nr ==1 && nnroi(1+nrois) >= 1  %map outside rois, this case never happern here.
                    vroi(i1,i2,i3,1) = 2;
                else % nr == 0   %do not map to the cortex
                    vroi(i1,i2,i3,1) = 1;
                end

            end
          end
      end
             
             roi_shtname = cell(1,nrois);
             for i=1:length(sroi_names)
                 roi_shtname{i} = strrm(sroi_names{i},['.' lr]);
                 roi_shtname{i} = strrm(roi_shtname{i},'.1D.roi');
             end
             if file_save
                 nSB = 3+(nrois+1)*3; %number of subbricks
                        
                 bl = 'Mapping property~roi index~Number of mapped rois~';
                 bl1 = sprintf('mapped to %s~',roi_shtname{:});
                 bl1o = 'mapped outside_rois~';
                 bl2 = sprintf('# of mapped nodes in %s~',roi_shtname{:});
                 bl2o = '# of mapped nodes outside_rois~'; 
                 bl3 = sprintf('is # of nodes maximum in %s~',roi_shtname{:});
                 bl3o = 'is # of nodes maximum outside_rois~';
                 
                 
                 info.BRICK_LABS = [bl,bl1,bl1o,bl2,bl2o,bl3,bl3o];
                 info.IDCODE_STRING = 'Malach';
       
                 info.DATASET_RANK(2) = nSB;  % the number of subbricks to save.
                 info.TYPESTRING = '3DIM_HEAD_FUNC';
                 info.SCENE_DATA = [0 0 1];  % first 0: orig view; second 0: 1 value; 1: matches TYPESTRING 
                 info.BRICK_STATS=[];  %minimum and maximum values of the subbrick;
                 info.BRICK_TYPES=3*ones(1,nSB);  %1: short, 3: float
                 info.BRICK_FLOAT_FACS = [];
           
                 info.IDCODE_DATE = date;
                 info.TAXIS_NUMS=[];
       
                 opt.Prefix = ['sroi2v2_',lr];
                 opt.OverWrite = 'y';
                 WriteBrik(vroi,info,opt);
                 disp(['Dataset ' opt.Prefix '+orig written to disk']); 
             end