function [mag,ph,max_mag,ph_mean,mag_4ph]=coil_combine(im,nois,consider_correlation,direct_add,do_scale,roi,calc_roi_only)
%  [mag,ph,max_mag,ph_mean]=coil_combine(im,nois,consider_correlation,direct_add,do_scale,roi)
% im: ro*pe*sl*nch*extra2
% nois: ro*nch
% assume no correlation between coils
% ph
% radian; 3/5/2017: changed to degree in int16 in units of 0.01 to save disk space
% mag: changed to uint16

if ~exist('calc_roi_only','var')
    calc_roi_only = false;
end

if ~exist('consider_correlation','var')
    consider_correlation=false;
end

if ~exist('direct_add','var')
   direct_add=false;  
end

if ~exist('do_scale','var')
    do_scale=true;
end

sz=size(im);
if length(sz)==4
   sz(5)=1; 
end

if sz(5)>1
    ph_mean=zeros(sz([1,2,3]),'int16');
end

ph=zeros(sz([1,2,3,5]),'int16');


mag=zeros(sz([1,2,3,5]),'single');

mag_4ph=zeros(sz([1,2,3,5]),'single');

if exist('nois','var')  && ~isempty(nois)
nois=nois-repmat(mean(nois,1),[size(nois,1),1]);
if ~consider_correlation
    tmp=diag(diag(nois'*nois));
    iRmat=inv(tmp);
else
    iRmat=inv(nois'*nois);

end

else
    iRmat=eye(size(im,4));
end

inois=diag(sqrt(iRmat));
%iRmat=eye(32);
roi=ri_d1(roi);
if direct_add

phc=angle(mean_roi(im(:,:,:,:,1),roi));
phc=reshape(phc,[1,1,1,length(phc)]);
phc=repmat(phc,[size(im,1),size(im,2),size(im,3),1]);
im_lowpass=abs(im(:,:,:,:,1)).*exp(1i*phc);

else
im_lowpass=lowPassHanningFilter2D(mean(im,5),64);
end            
                
for i=1:size(im,1)
    if mod(i,50)==0
       % disp(i);
    end
    for j=1:size(im,2)
        for k=1:size(im,3)
            
            if calc_roi_only && roi(i,j,k)==0
                continue;
            end
            var4pc=0;
            for l=1:sz(5)
                x=reshape(im(i,j,k,:,l),[sz(4),1]);
                %mag(i,j,k,iv)
                
                xref=reshape(im_lowpass(i,j,k,:),[sz(4),1]);
                
                imcmp=mean(conj(xref(:)).*x(:).*inois(:).^2);
                
                mag_4ph(i,j,k,l)=abs(sqrt(imcmp));
               
                imcmp2=mean(conj(x(:)).*x(:).*inois(:).^2);
                
                mag(i,j,k,l)=abs(sqrt(imcmp2));
                
                
                var4pc=var4pc+imcmp;% - make it consistent with the scanner;  2/15/17: 5:12 pm
                % 2/17/2017: + to match scanner; scanner sign changed because
                % setFreePa(1,2) was commented out.  Now they are uncommented.
                
                ph(i,j,k,l)=int16(angle(imcmp)/pi*18000);
                
            end
            if sz(5)>1
                ph_mean(i,j,k)=int16(angle(var4pc)/pi*18000);
            else
                ph_mean='';
            end
            
            
            
        end
        
        
    end
end

if do_scale
max_mag=max(mag(:));

mag=uint16(mag/max_mag*60000);

end