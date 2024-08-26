function WASSR_fit_varian(fname,mask,inc)
%WASSR_fit(fname,mask,inc)
if nargin==0
    help WASSR_fit
end

d=rdSdt(fname);

sz=size(d);

off=readPar(fname,'cfreq');
off=off(4:end);
if ~exist('inc','var')
    inc=1:length(off);
end

if ~exist('mask','var') || isempty(mask)
    mask=d(:,:,:,1)>max(d(:))*0.05;
end

d=reshape(d,[sz(1:3),length(off),sz(4)/length(off)]);
d=mean(d,5);
f=zeros(sz(1:3));
off2=off(inc)';

for i=1:size(d,1)
  disp(i);
    for j=1:size(d,2)
        for k=1:size(d,3)
            
            if mask(i,j,k)==0
                continue;
            end             
            y=d(i,j,k,inc);
            fitfunc=@(x) calcResidual(off2,squeeze(y),x);
            
            [fmin,ind]=min(y);
            if ind-1<5 || length(y)-ind<5
                continue;
            end
            
            finit=mean(off2(ind-1:ind+1));
            f(i,j,k)=fminsearch(fitfunc,finit);
        end
    end
end

writesdt4(f,[fname,'_centerF']);


function res=calcResidual(f2,y,x)


asym=z2asym_cutoff(f2,y,x);

res=sum(asym.^2)/length(asym);
