function [d,dim]=read_afni_sdt_images(varargin)

 % one_image: get only one sub-brik
 % brik: the brik to get if one_image = true;
 
 [d,dim]=ri(varargin{:});
 