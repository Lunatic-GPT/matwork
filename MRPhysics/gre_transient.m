function s=gre_transient(b,N,TR,T1,M0)
% s=gre_transient(b,N,TR,T1,M0)
% N is the number of TRs.
% b is the flip angle in degree; can be an array to alternate between flip
% angles
% M0: initial magnetization


if ~exist('M0','var')
    M0=1;
end
M=M0;
s=zeros(1,N);
for i=1:N
    
    fa=b(mod_1base(i,length(b)));
    
  s(i)=M*sin(fa*pi/180);
  M=relax(M*cos(fa*pi/180),TR,T1);
  
end

if nargout==0
    hold on;myPlot((1:N)*TR,s,'o-');
    ylabel('Signal Intensity (M_0)');
    xlabel('Time (s)');
end



function s=relax(m0,delay,T1)

s=1-(1-m0)*exp(-delay/T1);

