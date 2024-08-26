function k=fov_shift_kspace(k,dx,fov)
%k=fov_shift_kspace(k,dx,fov)

sz=size(k);
if numel(sz)==2 
    sz(3)=1;
end

n=size(k,1);
p1=exp(1i*2*pi*dx(1)/fov(1)*(1:n)');
p1=repmat(p1,[1,sz(2:end)]);

n=size(k,2);
p2=exp(1i*2*pi*dx(2)/fov(2)*(1:n));
p2=repmat(p2,[sz(1),1,sz(3:end)]);

k=k.*p1.*p2;