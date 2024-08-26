function [mag,ph,max_mag]=coil_combine_PC(im,nois,consider_correlation,do_scale)
%  [mag,ph]=coil_combine(im,nois,consider_correlation)
% im: ro*pe*sl*nch*extra1*extra2
% nois: ro*nch
% assume no correlation between coils
% ph: consider the phase difference between extra1=1 and 2:end; unit
% radian; 3/5/2017: changed to degree in int16 in units of 0.01 to save disk space
% mag: changed to uint16

if ~exist('consider_correlation','var')
    consider_correlation=false;
end

if ~exist('do_scale','var')
    do_scale=true;
end

sz=size(im);
if length(sz)==5
    sz(6)=1;
end
ph=zeros(sz([1,2,3,5,6])-[0,0,0,1,0],'single');


mag=zeros(sz([1,2,3,5,6]),'single');

if exist('nois','var') && ~isempty(nois)
    
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


%iRmat=eye(32);
tic;
for j=1:size(im,2)

    if mod(j,100)==0
        disp(j);
        toc;
    end
    
for k=1:size(im,3)
    
    for iv=1:sz(5)
        
        for l=1:sz(6)
         
            x=squeeze(im(:,j,k,:,iv,l));
            %mag(i,j,k,iv)
            %x=transpose(x(:));
            mag(:,j,k,iv,l)=diag(x*iRmat*x');
            
            if iv>1
                
                x1=squeeze(im(:,j,k,:,1,l));
                x2=squeeze(im(:,j,k,:,iv,l));
                
                
                ph(:,j,k,iv-1,l)=diag(x1*iRmat*x2');
                % - make it consistent with the scanner;  2/15/17: 5:12 pm
                % 2/17/2017: + to match scanner; scanner sign changed because
                % setFreePa(1,2) was commented out.  Now they are uncommented.
                
            
            end 
        end 
    end  
end
    
    if mod(j,100)==0
        tic;
    end
    
end

 mag=abs(sqrt(mag));
 ph=int16(angle(ph)/pi*18000);
 max_mag=max(mag(:));
 if do_scale


mag=uint16(mag/max_mag*60000);

 end
