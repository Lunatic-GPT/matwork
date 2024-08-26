function voxelSurfDist(f3dname,subject,thresh_ind,thresh,suma,sv,use_norm)
% voxelSurfDist(f3dname,subject,thresh_ind,thresh[,suma,sv,use_norm])
% This function finds the maximum, minimun and mean distances between nodes mapped to 
% the same voxel and append those distances as subbricks at the send of <f3dname> dataset.  
% 
%
% f3dname: file name of a 3d dataset with one sub-brick.
%          All the voxels with values greater than 0 will be considered
% subject: Name of the subject
% thresh_ind: the index of the subbrick used for thresholding voxels. 0
%             based.
% thresh: threshold value.  voxels with absolute value <= thresh
%         will not be considered.
% lr:      'lh' or 'rh'
% suma: suma directory. Default: '../../segmentation/<subject>/SUMA'
% sv:  surface volume aligned with the experiment. 
%      Default: '<subject>_SurfVol_Alnd_Exp+orig'
% use_norm: whether to use normals to Surf_A to generate the segments.
%      Default: true;
%
%      
%
% History: 2/05/2009, created by X.Z.
          
if ~exist('suma','var')
    suma = sprintf('../../segmentation/SUMA');
end

if ~exist('sv','var')
    sv = sprintf('%s_SurfVol_Alnd_Exp+orig',subject);
end

if ~exist('use_norm','var')
   use_norm = true;    
end

        
            [vdata,info]=BrikLoad(f3dname);
                        
            lr = 'lh';
            of1d = sprintf('temp_actvtn_v2s.%s.1D.dset',lr);
        
            if exist(of1d,'file') 
                delete(of1d);
            end
                        
            if ~use_norm
            cmd = sprintf(['3dVol2Surf' ... 
             ' -spec %s/%s_both.spec' ...
             ' -surf_A %s.smoothwm.asc'  ...
             ' -surf_B %s.pial.asc'  ...
             ' -sv %s' ...
             ' -grid_parent %s' ...  
             ' -map_func ave -f_steps 10' ... 
             ' -f_index voxels' ...
             ' -cmask ''-a %s[%d] -expr step(a-%f)'' ' ...
             ' -out_1D %s'], ...
              suma,subject,lr,lr,sv,f3dname,f3dname,thresh_ind,thresh,of1d);
            else 
             cmd = sprintf(['3dVol2Surf_local' ... 
             ' -spec %s/%s_both.spec' ...
             ' -surf_A %s.smoothwm.asc'  ...
             ' -sv %s' ...
             ' -use_norms' ...
             ' -norm_len 3' ...
             ' -grid_parent %s' ...  
             ' -map_func max -f_steps 10' ... 
             ' -f_index voxels' ...
             ' -cmask ''-a %s[%d] -expr step(a-%f)'' ' ...
             ' -out_1D %s'], ...
              suma,subject,lr,sv,f3dname,f3dname,thresh_ind,thresh,of1d);
            end
            
             disp('Converting from surface to volume');
             unix(cmd);
             disp(sprintf('File %s written to disk',of1d));
             
             sdata = textread(of1d,'','commentstyle','shell');            
             
             % find the nodes that are mapped from each voxel
             
             
            % loop over the active voxels
             nodes = cell(size(vdata,1),size(vdata,2),size(vdata,3));  % store nodes mapped to each voxel
             
             data = zeros(size(vdata,1),size(vdata,2),size(vdata,3),6);
         
             for row = 1:size(sdata,1)
                 i = sdata(row,3) + 1;
                 j = sdata(row,4) + 1;
                 k = sdata(row,5) + 1;
                nodes{i,j,k}(end+1) =sdata(row,1); 
             end
             
             for i=1:size(data,1)
                 for j=1:size(data,2)
                    for k=1:size(data,3)
                        data(i,j,k,1) = length(nodes{i,j,k});         % number of nodes.
                    end
                 end
              end
            disp('calculate the distances between nodes mapped to each voxel');
             avind = find(vdata(:,:,:,thresh_ind+1)>thresh & data(:,:,:,1)>1);  % find active voxels
             
             debug_avind = find(vdata(:,:,:,thresh_ind+1)<thresh & data(:,:,:,1)>0,1);
             if ~isempty(debug_avind)
                 disp('You should not see this.  Check your program');
             end
             
             nav = length(avind);
             for i=1:nav
                [a,b,c]= ind2sub(size(vdata),avind(i));
                nn = data(a,b,c,1);
                pairs = zeros(nn*(nn-1)/2,2); % node pairs
                ipar = 0; % pair index
                
                %generate the node pairs 
                for j = 1:nn-1
                    for k = j+1:nn
                        ipar = ipar+1;
                        pairs(ipar,1) =nodes{a,b,c}(j);
                        pairs(ipar,2) =nodes{a,b,c}(k);
                        
                    end
                end
                
                fid = fopen('temp_node_pairs.1D','w');
                
                for ip = 1:ipar
                fprintf(fid,'%d\t%d\n',pairs(ip,1),pairs(ip,2));
                end
                fclose(fid);
                
                % call SurfDist
                cmd = sprintf('SurfDist -i %s/%s.smoothwm.asc -input temp_node_pairs.1D > distance.txt',suma,lr);
                               
                unix(cmd); 
                dist = textread('distance.txt','','commentstyle','shell');
                data(a,b,c,2) = max(dist(:,3));
                data(a,b,c,3) = min(dist(:,3));
                data(a,b,c,4) = mean(dist(:,3));
                          
                           
                if(mod(i,5)==0)
                    disp(sprintf('%d of %d voxels processed',i,nav));
                end
                 
             end
                 maxdata = data(:,:,:,2);
                 meandata = data(:,:,:,4);
                 maxmax = max(maxdata(:));
                 maxmean = max(meandata(:));
                 
                 for i=1:size(data,1)
                     for j=1:size(data,2)
                         for k=1:size(data,3)
                             if data(i,j,k,2) > 0;
                                data(i,j,k,5) = maxmax-data(i,j,k,2);
                             end
                             if data(i,j,k,4) > 0;
                                data(i,j,k,6) = maxmean-data(i,j,k,4);
                             end
                         end
                     end
                 end
                 
                            
                 nSB = size(vdata,4)+6;
                 brikdata = zeros(size(vdata,1),size(vdata,2),size(vdata,3),nSB);
                 brikdata(:,:,:,1:size(vdata,4)) = vdata;
                 brikdata(:,:,:,size(vdata,4)+1:nSB) = data;
                 
                 bl = '~number of nodes~max distance~min distance~mean distance~inverse max~inverse mean~';
                 
                 info.BRICK_LABS = [info.BRICK_LABS,bl];
                 info.IDCODE_STRING = [info.IDCODE_STRING '&Vol to Surf mapping'];
       
                 info.DATASET_RANK(2) = nSB;  % the number of subbricks to save.
                 info.TYPESTRING = '3DIM_HEAD_FUNC';
                 info.SCENE_DATA = [0 0 1];  % first 0: orig view; second 0: 1 value; 1: matches TYPESTRING 
                 info.BRICK_STATS=[];  %minimum and maximum values of the subbrick;
                 info.BRICK_TYPES=3*ones(1,nSB);  %1: short, 3: float
                 info.BRICK_FLOAT_FACS = [];
           
                 info.IDCODE_DATE = date;
                 info.TAXIS_NUMS=[];
       
                 %opt.Prefix = 'retino_map_phase2+orig';
                 %opt.Prefix = out_prefix;
                 out_prefix = strrm(f3dname,'+orig');
                 opt.Prefix = [out_prefix,'_SurfDist'];
                 opt.OverWrite = 'y';
                 WriteBrik(brikdata,info,opt);
                 disp(['Dataset ' opt.Prefix '+orig written to disk']); 
            