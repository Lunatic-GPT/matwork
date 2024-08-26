function fx = dct_trunc(x,nseg)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if mod(length(x),nseg)~=0
    error('mod(length(x),nseg)~=0');
end

fx=zeros(size(x));
l=length(x)/nseg;

  for i=1:nseg
    fx((i-1)*l+1:i*l)=dct(x((i-1)*l+1:i*l));
  end

end

