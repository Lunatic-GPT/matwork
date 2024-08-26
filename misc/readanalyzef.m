function [d,pixdim,dtype] = readanalyzef(filename,vol)

for i=1:length(vol)
    [d(:,:,:,i),pixdim]=readanalyze(filename,vol(i));
end
