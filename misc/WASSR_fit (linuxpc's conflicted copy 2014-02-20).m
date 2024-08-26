function WASSR_fit(fname,mask,inc,frange,off)
%WASSR_fit(fname[,mask,inc,frange,off])
% inc: index to include in the fit, 1 based.
% frange: +- frequency range from minimum to include in the fit. Units: Hz.
% set inc to [] for this option.
if nargin==0
    help WASSR_fit
end

if exist('frange','var')
    inc=[];
end

if exist([fname,'+orig.HEAD'],'file')
[d,info]=BrikLoad([fname,'+orig']);
elseif exist([fname,'_recon_phasecorr+orig.HEAD'],'file')
    [d,info]=BrikLoad([fname,'_recon_phasecorr+orig']);
else
    
    [d,info]=BrikLoad([fname,'_recon+orig']);
end    

sz=size(d);

if ~exist('off','var')
  off=readbPar([fname,'/method'],'freqlist');
end

if ~exist('inc','var')  || isempty(inc)
    inc=1:length(off);
end

if ~exist('mask','var') || isempty(mask)
    mask=(d(:,:,:,1)>max(d(:))*0.02)|(d(:,:,:,end)>max(d(:))*0.02);
elseif ~isempty(strfind(mask,'.mat'))
    mask=load(mask);
    mask=mask.roi;
else
    mask=BrikLoad(mask);
    
end

d=reshape(d,[sz(1:3),length(off),sz(4)/length(off)]);
d=mean(d,5);
f=zeros(sz(1:3));
write_afni(d,[fname,'_zspectrum'],info);

for i=1:size(d,1)
  fprintf('%d ',i);
    for j=1:size(d,2)
        for k=1:size(d,3)
       %     disp([i,j,k]);
            if mask(i,j,k)==0
                continue;
            end
            
             
            if i==20 && j==31 
  %              disp('');
            end
            if exist('frange','var')
            
                y=squeeze(d(i,j,k,:));
                
                [tmp,ind]=min(y);
                
                ind2=find(abs(off(ind)-off)<=frange);
               f(i,j,k)=WASSR_fit_1d(off(ind2),y(ind2));

            else
             y=squeeze(d(i,j,k,inc));
             off2=off(inc);
            
              f(i,j,k)=WASSR_fit_1d(off2,squeeze(y));
            end
        end
    end
end
disp('');
write_afni(f,[fname,'_centerFnew'],info);


function res=calcResidual(f2,y,x)


asym=z2asym_cutoff(f2,y,x);

res=sum(asym.^2)/length(asym);
