function afni2sdt(prefix,keep_file)
%afni2sdt(prefix,keep_file)
% first dim: ro; second dim: pe; third dim: slice;
% 1/22/2011: Wrote it. Tested for single-slice 2D data. XZ
% nophase: do not save phase image.default true.

if ~exist('keep_file','var')
    keep_file=true;
end

a=BrikLoad([prefix,'+orig']);

writesdt4(a,prefix);

if ~keep_file
    
    delete([prefix,'+orig.*']);
end



