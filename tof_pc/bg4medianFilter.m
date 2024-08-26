function [bg_mean,bg_std,bg]=bg4medianFilter(dname,npix,m_bg)
% background is estimated from the mean along the time points after median filtering (4th
% dimension)
%



a=ri(dname);

a=double(a);
bg=a(:,:,:,1)*0;

%for nt=1:size(a,4)
    for sl=1:size(a,3)
        
        a_tmp=mean(a(:,:,sl,:),4);      
        bg(:,:,sl)=median_filter(a_tmp,npix,m_bg(:,:,sl));
             
    end
%end

%dnew=a-repmat(bg,[1,1,1,size(a,4)]);

    bg_mean=mean(bg(m_bg>0));
    bg_std=std(bg(m_bg>0));

function res=median_filter(d,l,m_wm)

res=d.*0;

if mod(l,2)~=1
    l=l+1;
  %  warning('Use an odd size for the median filter: size = %d',l);
end

sz=size(d);
for i=(l+1)/2:sz(1)-(l-1)/2
    
    for j=(l+1)/2:sz(2)-(l-1)/2
        if m_wm(i,j)==0
            continue;
        end
        
        indi=i-(l-1)/2:i+(l-1)/2;
        
        indj=j-(l-1)/2:j+(l-1)/2;
        
        
        x=d(indi,indj);
        sel=m_wm(indi,indj);
        res(i,j)=median(x(sel>0));
        
    end
end


