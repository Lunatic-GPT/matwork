function res = Wavelet_ISU(N12,level)
% res = Wavelet_ISU(N12,level)
%
% implements a wavelet operator
%N12 is the size of the 2D/3D images.
if ~exist('level','var')
    level=2;
end
dwtmode('per');

N=N12(1);
 Xtmp=eye(N);
 W=zeros(N);
  for p=1:N
      W(:,p)=wavedec(Xtmp(:,p),level,'db4');
  end
res.wavl = W;

if numel(N12)==2 && N12(1)~=N12(2) 
 N=N12(2);
 Xtmp=eye(N);
 W2=zeros(N);
  for p=1:N
      W2(:,p)=wavedec(Xtmp(:,p),level,'db4');
  end
 res.wavr = W2;
else
 res.wavr = W;
end

if numel(N12)==3 && N12(1)~=N12(3) 
 N=N12(3);
 Xtmp=eye(N);
 W3=zeros(N);
  for p=1:N
      W3(:,p)=wavedec(Xtmp(:,p),level,'db4');
  end

  res.wavs = W3;
else
  res.wavs = W;
end

res.N=N12;

res.adjoint = 0;
%res.qmf = MakeONFilter(filterType, filterSize);
res = class(res,'Wavelet_ISU');
