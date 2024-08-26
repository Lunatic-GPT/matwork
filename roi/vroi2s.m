function sroi = vroi2s(f3dname,mask_ind,out_prefix,subject,lr,suma,sv)
% 
% sroi = vroi2s(f3dname,mask_ind,out_prefix,subject,lr[,suma,sv])
% This function maps rois defined on the surface to 3d volume.           
% Because of the non-unique nature of the surface to volume nature. 
% voxels can be classified into the following three types:
%
%     2. nodes that map to a single ROI
%     1. nodes that map to multiple ROIs.
%     0. nodes that do not map to any ROIs.
%     
% f3dname: file name of a 3d dataset with one sub-brick.
%          All the voxels with values greater than 0 will be considered
% mask_ind: Index of the mask brick.  1 based.
% out_prefix: prefix of the output file <out_prefix>.<lr>.1D.dset
% subject: Name of the subject
% lr:      'lh' or 'rh'
% suma: suma directory. Default: '../../../segmentation/<subject>/SUMA'
% sv:  surface volume aligned with the experiment. 
%      Default: '<subject>_SurfVol_Alnd_Exp+orig'
%
% sroi is a matrix with the following columns 
% 1. Node index
% 2. mapping property of the node with values of 0-2 as specified above
% 3. number of mapped rois.   
% 4. number of mapped roi voxels.
% 5. node values. One column for each subbrick.
%    Maximum value of all mapped voxels.
%
% Note: This function only process one hemisphere. 
%       
%
% History: 1/29/2009, created by X.Z.
          
if ~exist('suma','var')
    suma = sprintf('../../../segmentation/%s/SUMA',subject);
end

if ~exist('sv','var')
    sv = sprintf('%s_SurfVol_Alnd_Exp+orig',subject);
end


       
            [vdata,info]=BrikLoad(f3dname);
                        
            nf = size(vdata,1);
            np = size(vdata,2);
            ns = size(vdata,3);
            nbricks = size(vdata,4);
                                          
            of1d = sprintf('temp_actvtn_v2s.%s.1D.dset',lr);
            of1d_seg = sprintf('temp_segcoord.%s.1D.dset',lr); 
            
            if exist(of1d,'file') 
                delete(of1d);
            end
            if exist(of1d_seg,'file') 
                delete(of1d_seg);
            end
             
            cmd = sprintf(['3dVol2Surf' ... 
             ' -spec %s/%s_both.spec' ...
             ' -surf_A %s.smoothwm.asc'  ...
             ' -surf_B %s.pial.asc'  ...
             ' -sv %s' ...
             ' -grid_parent %s' ...  
             ' -map_func max -f_steps 10' ... 
             ' -f_index voxels' ...
             ' -save_seg_coords %s' ...
             ' -cmask ''-a %s[%d] -expr step(a-0.01)'' ' ...
             ' -out_1D %s'], ...
              suma,subject,lr,lr,sv,f3dname,of1d_seg,f3dname,mask_ind-1,of1d);
             disp('Coverting from surface to volume');
             disp(sprintf('File %s written to disk',of1d));
             unix(cmd);
             
             sdata = textread(of1d,'','commentstyle','shell');
             nc = textread(of1d_seg,'','commentstyle','shell');
             
             % find the nodes that are mapped from each voxel
             
             nnodes = size(sdata,1);
             sroi = zeros(nnodes,4+nbricks);
            % loop over the surface nodes.
             for i=1:nnodes
                
                nr = 0;  % belongs to how many rois.
                indr = []; % index of rois.
                
                for ivox = 1:sdata(i,6)
                                              
                    ijk = xyz2ijk(nc(i,2+3*(ivox-1):1+3*ivox),info);
                    if isempty(indr) || ...
                       isempty(find(vdata(ijk(1),ijk(2),ijk(3))==indr,1))
                       nr = nr+1;
                       indr(end+1) = vdata(ijk(1),ijk(2),ijk(3));
                    end
                    
                end
 
                sroi(i,1) = sdata(i,1); %node index
                sroi(i,4) = sdata(i,6); %number of voxels
                sroi(i,3) = nr;
                sroi(i,5:4+nbricks) = sdata(i,7:6+nbricks);
                
                if nr == 1  %map to a single roi
                    sroi(i,2) = 2;
                elseif nr > 1  % map to more than one roi.
                    sroi(i,2) = 1;
                else % do not map to roi. This should not occur.
                    disp('Program error');
                    sroi(i,2) = 0;
                end
                
             end
             
             fname = sprintf('%s.%s.1D.dset',out_prefix,lr);
             fid = fopen(fname,'w');
             fprintf(fid,'#node index\tnode type\tnumber of rois\tnumber of voxels\t');
             fprintf(fid,'%s\n',info.BRICK_LABS);                
             for i=1:size(sroi,1)
                
                fprintf(fid,'%d\t', sroi(i,1:4));
                fprintf(fid,'%3.2f\t',sroi(i,5:end));
                fprintf(fid,'\n');
             end
             fclose(fid);
             disp(['FILE: ' fname, ' written to disk']);

function ijk = xyz2ijk(xyz,info)
% ijk indices in afni.  
% This is one less than the indices in matlab.    
                dxyz = info.DELTA;
                orgn = info.ORIGIN; 
                ijk = floor((xyz-orgn)./dxyz+0.4999);