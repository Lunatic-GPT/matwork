function [cv_all,cv_all_same,cv_all_nobg,cv_all_bg]=retro_time_course(phase_file,pvs_mask,wm_mask,venc,use_peak,neighborWidth,phase2deg,mdfSize,vsign)
%[v_all,v_all_same,v_all_nobg,v_all_bg]=retro_time_course(phase_file,pvs_mask,wm_mask,venc,use_peak,neighborWidth,phase2deg[,mdfSize])
% file: phase in units of degree
% if mag file: venc should be 'tof', otherwise the venc value
% mdfSize: if mdfSize is given, median filter will be performed on the
% phase data to estimate the background.
if ~exist('phase2deg','var')
    phase2deg=1;
end

if isa(pvs_mask,'char')
    try
    mall=ri(pvs_mask,[],[],'d');
    catch
        mall=ri(pvs_mask);
    end
else
    mall=pvs_mask;
end

if ~exist('vsign','var')
    vsign=ones(1,size(mall,4));
end

if isa(phase_file,'char')
    try
      b=ri(phase_file);
    catch
      b=ri(phase_file);
    end
else
    b=phase_file;
end

if isempty(wm_mask)
    m_wm=ones(size(mall));
elseif isa(wm_mask,'char')
    try
      m_wm=ri(wm_mask,[],[],'d');
    catch
      m_wm=ri(wm_mask);
    end
else
    m_wm=wm_mask;
end

if ~exist('do_plot','var')
do_plot=false;
end

    b=double(b)*venc/180*phase2deg;

    vsign=reshape(vsign,[1,1,length(vsign)]);
    b=b.*repmat(vsign,[size(b,1),size(b,2),1,size(b,4)]);
    
%use_peak=true;  %use only the voxel with the maximum phase in each roi
%% assume different baseline


for is=1:size(mall,3)
v_all=[];    % assume different baseline
v_all_same=[];  % same baseline
v_all_nobg=[];  % no baseline subtraction

v_all_bg=[];
%ncol=ceil(sqrt(max(m(:))));

%nrow=ceil(max(m(:))/ncol);

m=clusterize2(mall(:,:,is),1);

for i=1:max(m(:))

  %  subplot(nrow,ncol,i);
    
 %   m_bg=roiRing(m==i,ringWidth);
    m_bg=roiCOMNeighbor(m==i,neighborWidth);
    
    m_bg=m_bg&m_wm(:,:,is)&m==0;
    
   
    
   fprintf('Background ROI voxels = %d\n',sum(m_bg(:)));
    
   
   if ~exist('mdfSize','var')  || isempty(mdfSize)
       v_bg=mean_roi(b(:,:,is,:),m_bg);
      
   else
       v_bg=bg4medianFilter(b(:,:,is,:),mdfSize,m_bg);
   end
   
           %   disp(v_bg);
    if use_peak
    
        mb=mean(b(:,:,is,:),4);
        maxb=max(mb(m==i));
     %   ind=find();
        disp(i);
        mpeak=m*0;
        mpeak(m==i&mb==maxb)=1;
       
        v=mean_roi(b(:,:,is,:),mpeak)-v_bg;
        v_same=mean_roi(b(:,:,is,:),mpeak)-mean(v_bg);
        fprintf('v_bg mean = %f\n',mean(v_bg));
        v_nobg=mean_roi(b(:,:,is,:),mpeak);
        
    else
           v=mean_roi(b(:,:,is,:),m==i)-v_bg;
          v_same=mean_roi(b(:,:,is,:),m==i)-mean(v_bg);
          v_nobg=mean_roi(b(:,:,is,:),m==i); 
       
    
    end
    
    v_all_same=cat(1,v_all_same,v_same(:)');
    v_all=cat(1,v_all,v(:)');
    v_all_nobg=cat(1,v_all_nobg,v_nobg(:)');
    v_all_bg=cat(1,v_all_bg,v_bg(:)');
    
end

if size(mall,3)==1
    cv_all_same=v_all_same;
    cv_all=v_all;
    cv_all_nobg=v_all_nobg;
    cv_all_bg=v_all_bg;
else
       cv_all_same{is}=v_all_same;
    cv_all{is}=v_all;
    cv_all_nobg{is}=v_all_nobg;
    cv_all_bg{is}=v_all_bg;
    
end

end

