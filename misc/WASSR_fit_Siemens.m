function WASSR_fit_Siemens(fname,mask,off,inc)
%WASSR_fit(fname,mask,off,inc)
% inc: index to include in the fit, 1 based.
% off: frequency list. Units: ppm.
%mask: can be [], then will be calculated based on 0.02*max(d(:))

if nargin==0
    help WASSR_fit_Siemens
end


d=ri(fname);

sz=size(d);


if ~exist('inc','var')  || isempty(inc)
    inc=1:length(off);
end

if ~exist('mask','var') || isempty(mask)
    mask=(d(:,:,:,1)>max(d(:))*0.02)|(d(:,:,:,end)>max(d(:))*0.02);
else ~isempty(strfind(mask,'.mat'))
    mask=ri(mask);
    
end

d=reshape(d,[sz(1:3),length(off),sz(4)/length(off)]);
d=mean(d,5);
f=zeros(sz(1:3));

for i=1:size(d,1)
fprintf('%d ',i);

if mod(i,20)==0
    disp(' ');
end
ftmp=zeros(sz(2:3));
    for j=1:size(d,2)
        for k=1:size(d,3)
           
            if mask(i,j,k)==0
                continue;
            end
            
             y=squeeze(d(i,j,k,inc));
             off2=off(inc);
            
              ftmp(j,k)=WASSR_fit_1dt(off2(:),squeeze(y));
       end
       
    end
     f(i,:,:)=ftmp;
end
disp(' ');
fname=strtok(fname,'.');
save([fname,'_B0'],'f');


function res=calcResidual(f2,y,x)


asym=z2asym_cutoff(f2,y,x);

res=sum(asym.^2)/length(asym);
