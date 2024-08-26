function fl_fq_retro_medianFilter(dname,npix,m_wm)
% fl_fq_pc(dname,npix,m_wm)
%

    prefix=dname;
    dir_str=dir(fullfile(dname,'*_P_*'));
if length(dir_str)==1
    dname=fullfile(dname,dir_str(1).name);
end


a=ri(dname,1);

a=single(a);

a_max=max(a(:));
a_min=min(a(:));

a=(a-(a_min+a_max)/2)/(a_max-a_min)*360;

a_samebg=a*0;

a_indvbg=a*0;
bg_mean=zeros(size(a(:,:,:,1)));

bg=a*0;

if ~exist('m_wm','var')
    m_wm=ones(size(a(:,:,:,1)));
else
    try 
    m_wm=ri(m_wm,'','','d');
    catch
        m_wm=ri(m_wm);
    end
end

for sl=1:size(a,3)
    
    a_tmp=mean(a(:,:,sl,:),4);
    bg_mean(:,:,sl)=median_filter(a_tmp,npix,m_wm);
    
    a_samebg(:,:,sl,:)=(a(:,:,sl,:)-repmat(bg_mean(:,:,sl),[1,1,1,size(a,4)])).*repmat(m_wm(:,:,sl),[1,1,1,size(a,4)]);
    
    if size(a,4)>1
    for k=1:size(a,4)
    
        disp([sl,k]);
        a_tmp=mean(a(:,:,sl,k),4);
        bg(:,:,sl,k)=median_filter(a_tmp,npix,m_wm);
        
        a_indvbg(:,:,sl,k)=(a(:,:,sl,k)-bg(:,:,sl,k)).*m_wm(:,:,sl);
        
    end
    end
    
end
%%

if length(prefix)>4 && strcmp(prefix(end-3:end),'.mat')
    prefix=prefix(1:end-4);
else
    prefix=strtok(prefix,'.');
end

if size(a_indvbg,4)>1
    
    save([prefix,'_indvbg_mft.mat'],'a_indvbg');
    save([prefix,'_samebg_mft.mat'],'a_samebg');
    ma=mean(a_samebg,4);
    save([prefix,'_samebg_mft_mean.mat'],'ma');
    save([prefix,'_bg.mat'],'bg');
    save([prefix,'_bg_mean.mat'],'bg_mean');
     
else
    save([prefix,'_mft.mat'],'a_samebg');
    
    save([prefix,'_bg.mat'],'bg_mean');
    
end


function res=median_filter(d,l,m_wm)

res=d.*0;

if mod(l,2)~=1
    error('Please use an odd size for the median filter');
end

sz=size(d);
for i=(l+1)/2:sz(1)-(l-1)/2
   disp(i);
    for j=(l+1)/2:sz(2)-(l-1)/2
        indi=i-(l-1)/2:i+(l-1)/2;
        
        indj=j-(l-1)/2:j+(l-1)/2;
        %{
        indi=i-(l-1)/2:i+(l-1)/2;
        indi(indi<1)=indi(indi<1)+size(d,1);
        indi(indi>size(d,1))=indi(indi>size(d,1))-size(d,1);
        
        indj=j-(sz-1)/2:j+(sz-1)/2;
        indj(indj<1)=indj(indj<1)+size(d,2);
        indj(indj>size(d,2))=indj(indj>size(d,2))-size(d,2);
        %}
    
        
        if m_wm(i,j)==0
            continue;
        end
        x=d(indi,indj);
        sel=m_wm(indi,indj);
        res(i,j)=median(x(sel>0));
        
    end
end


