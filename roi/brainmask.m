function m=brainmask(d,threshold,nvox)
% m=brainmask(d,threshold,nvox)
 
 
if ~exist('nvox','var')
    nvox=200;
end


m=d>threshold;

m=clusterize2_2d(m,nvox,2);

for i=1:size(m,3)
  m(:,:,i)=imfill(m(:,:,i),'holes');
end