function [nv,val]=nv_roi(data,mask)
% [nv,val]=nv_roi(data,mask)
% do average over roi defined in mask.
% data can be file name or 4d matrix
% mask can be file name or 3d matrix

tic;

if isa(data,'char') && isa(mask,'char')
    cmd = sprintf('3dmaskave -mask ''%s'' ''%s'' > maskave_temp.1D',mask,data);
    unix(cmd);
    
    [ts_tmp,nvox_temp,temp_str] = textread('maskave_temp.1D','%f [%d %s');
    if ~strcmp(temp_str,'voxels]')
        error('3dmaskave output file error');
    end
    if numel(nvox_temp) == 0
        nv = 0;
        val = 0;
    else
        nv = nvox_temp(1);
        val = ts_tmp;
    end
    
else
    if isa(data,'char')
        data = BrikLoadf(data);
    elseif isa(mask,'char')
        mask = BrikLoad(mask);
    end
    data=squeeze(data);
    mask=squeeze(mask);
      if ndims(data)- ndims(mask)~=1   && ndims(data)- ndims(mask)~=0 
        error('ndims(data)- ndims(mask)~=1  && ndims(data)- ndims(mask)~=0');
      end
    
    if ndims(data)==3 && ndims(mask)==2
        data =reshape(data,[size(data,1),size(data,2),1,size(data,3)]);
        mask =reshape(mask,[size(mask,1),size(mask,2),1,size(mask,3)]);
    end
    
    val = zeros(1,size(data,4));
    for i=1:size(data,4)
      tmp = data(:,:,:,i);
      val(i) = mean(tmp(mask>0));
    end
    
    nv = length(find(mask>0));

end

disp([mfilename ' finish in ', num2str(toc), ' s']);