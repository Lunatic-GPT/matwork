function d=setv_ind(d,ind,val)
% ind: 1D index of d(:,:,:,1);
% val: vec of size 1*size(d,4)

sz=size(d);

d=reshape(d,[sz(1)*sz(2)*sz(3),sz(4)]);



for i=1:length(val)
  
d(ind,i)=val(i);
end

d=reshape(d,sz);
