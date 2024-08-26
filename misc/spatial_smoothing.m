function a_new = spatial_smoothing(a,sig,prefix)
%a_new = spatial_smoothing(a,sig,prefix)
%sig: 1*3 matrix in units of voxel size.
% support 4D data.  If prefix not given, then do not save.
% FWHM=sig*2*sqrt(2*log(2));

if isa(a,'char')
    prefix=strtok(a,'.');
    prefix=[prefix,'_smooth'];
    a=ri(a);
end

n=ceil(sig*5)*2;


sig(sig==0)=1;
n(n==0)=1;

if size(a,3)==1
[x,y]=meshgrid(-n(1)/2+1:n(1)/2,-n(2)/2+1:n(2)/2);
g=exp(-(x.^2/sig(1)^2/2+y.^2/sig(2)^2/2));
    
else
    
[x,y,z]=meshgrid(-n(1)/2+1:n(1)/2,-n(2)/2+1:n(2)/2,-n(3)/2+1:n(3)/2);
g=exp(-(x.^2/sig(1)^2/2+y.^2/sig(2)^2/2+z.^2/sig(3)^2)/2);

end

g=g/sum(g(:));
sz=size(a);
if ndims(a)==2
    sz(3)=1;
end
a_new=zeros(sz);
for i=1:size(a,4)
    for j=1:size(a,5)
        a_new(:,:,:,i,j)=convn(a(:,:,:,i,j),g,'same')/sum(g(:));
       % a_new(:,:,:,i,j)=b(ceil(n(1)/2):ceil(n(1)/2)+sz(1)-1,ceil(n(2)/2):ceil(n(2)/2)+sz(2)-1,ceil(n(3)/2):ceil(n(3)/2)+sz(3)-1);
    end
end


if exist('prefix','var')
    save([prefix,'.mat'],'a_new');
end


