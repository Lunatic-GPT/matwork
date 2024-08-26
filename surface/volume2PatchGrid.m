function zi = volume2PatchGrid(f3dname,subject,patch_name,thresh_ind,thresh,lr,suma,sv)
% zi = volume2PatchGrid(f3dname,subject,patch_name,lr,[,suma,sv])
% converts data from 3d to an evenly spaced grid centered around a flat patch. 
%
% f3dname: file name of a 3d dataset.
%          All the voxels with values greater than 0 will be considered
% subject: Name of the subject
% patchname: the middle part of the patch file.  
%       The file name should be [r|l]h.<patchname>  
% thresh_ind: the index of the subbrik to be used as threshold. 0 based
% thresh:   threshold
% lr: hemisphere, lh or rh.
% suma: suma directory. Default: '../../../segmentation/<subject>/SUMA'
% sv:  surface volume aligned with the experiment. 
%      Default: '<subject>_SurfVol_Alnd_Exp+orig'
%
%      
%
% History: 2/14/2009, created by X.Z.
          imSize = 256;
          irange = 1;  % do interpratation for the current pixel 
                       % if there is a non-zero voxel within irange pixels
                       % from the current one.
if ~exist('suma','var')
    suma = sprintf('../../segmentation/SUMA');
end

if ~exist('sv','var')
    sv = sprintf('%s_SurfVol_Alnd_Exp+orig',subject);
end


            of1d = sprintf('temp_actvtn_v2s.%s.1D.dset',lr);
        
            if exist(of1d,'file') 
                delete(of1d);
            end
                        
            cmd = sprintf(['3dVol2Surf' ... 
             ' -spec %s/%s_both.spec' ...
             ' -surf_A %s.smoothwm.asc'  ...
             ' -sv %s' ...
             ' -use_norms' ...
             ' -norm_len 3' ...
             ' -grid_parent %s' ...  
             ' -map_func ave -f_steps 10' ... 
             ' -f_index voxels' ...
             ' -skip_col_1dindex' ...  %    : do not output 1dindex column
             ' -skip_col_i' ... %            : do not output i column
             ' -skip_col_j' ... %            : do not output j column
             ' -skip_col_k' ... %            : do not output k column
             ' -skip_col_vals' ...    %     : do not output vals column
             ' -cmask ''-a %s[%d] -expr step(a-%f)''' ... 
             ' -out_1D %s'], ...
              suma,subject,lr,sv,f3dname,f3dname,thresh_ind,thresh,of1d);
           
             disp('Converting from surface to volume');
             unix(cmd);
             disp(sprintf('File %s written to disk',of1d));
             
             sdata = textread(of1d,'','commentstyle','shell');            
             
             patchf = sprintf('%s/%s.%s',suma,lr,patch_name); 
             [nodes, coords] = readPatch(patchf);
             
                       
               x = zeros(size(sdata,1),1);
               y = zeros(size(sdata,1),1);
               z = zeros(size(sdata,1),size(sdata,2)-1);
               nnzero = 0;
     % read in the data.
               for in = 1:size(nodes,2)
                   
                   row = find(sdata(:,1)==nodes(in),1);
                   if ~isempty(row)
                       nnzero = nnzero+1;
                      x(nnzero) = coords(1,in);
                      y(nnzero) = coords(2,in);
                      z(nnzero,1:end) = sdata(row,2:end); 
                   end
               end
              
               xmin = min(coords(1,:));
               xmax = max(coords(1,:));
               ymin = min(coords(2,:));
               ymax = max(coords(2,:));
               zmin = min(z,[],1);
               zmax = max(z,[],1);
               
               xstep = (xmax-xmin)/(imSize-1);
               ystep = (ymax-ymin)/(imSize-1);
               xi = xmin:xstep:xmax;
               yi = ymin:ystep:ymax;
               zi = zeros(imSize,imSize,size(z,2));
               
             % set value =-1 for pixels outside the activated region 
               nzi = zeros(imSize,imSize);
               xstep2 = (xmax-xmin)/imSize;
               ystep2 = (ymax-ymin)/imSize;
             
              for i=1:nnzero
                
                  xind = floor((x(i)-xmin)/xstep2)+1;
                  yind = floor((y(i)-ymin)/ystep2)+1;
                  
                  if xind > imSize 
                      xind = imSize;
                  end
                  if yind > imSize
                      yind = imSize;
                  end
                  
                 if z(i,thresh_ind+1)>0 
                   nzi(yind,xind) = nzi(yind,xind)+1;
                 end
              end
            
              for i=1:imSize
                  for j=1:imSize
                      if nzi(i,j)==0
                          
                         if isempty(find( (xi(j) - x(1:nnzero)).^2 + (yi(i)-y(1:nnzero)).^2<(irange*xstep2)^2,1) )
                             nnzero = nnzero+1;
                             x(nnzero) = xi(j);
                             y(nnzero) = yi(i);
                             for k=1:size(z,2)
                              z(nnzero,k)=NaN; %-(zmax(k)-zmin(k))*0.01;
                             end
                         end
                      end
                      
                  end
              end
                
        
      for i = 1:size(z,2)
        zi(:,:,i) = griddata(x(1:nnzero),y(1:nnzero),z(1:nnzero,i),xi,yi');
      end
               
      % one can use the following to set the un-activated regions look
      % gray.
      
      % colordata = colormap;
      % colordata(1,:) = [0.5,0.5,0.5];
      % colormap(colordata);
      
      
               
                   
               
                
       
            