function fl_fq_retro_pc(dname,mask_factor,roi_name,inv_off)
% fl_fq_pc(dname,mask_factor,mag_name)
% mask_factor: intensity threshold scale; threshold =
% mask_factor*max_intensity; default: 0.05.
% roi_name: can be a mask file/matrix or a dicom folder.  In the later, mask will
% be generated from 0.1*max(d(:)).  
% If not given, then assume the phase file is under dname/*_P_* and roi_name is under dname
% results in degrees
if ~exist('mask_factor','var')
    mask_factor=0.1;
end

if ~exist('inv_off','var')
    inv_off = false;
end

if ~exist('roi_name','var') || isempty(roi_name)

    roi_name=dname;
    dir_str=dir(fullfile(dname,'*_P_*'));
    dname=fullfile(dname,dir_str(1).name);
    
end

prefix=dname;
if isa(roi_name,'char')
    
    im=ri(roi_name,1);
    im2=mean(im,4);
    m=clusterize2(im2>mask_factor*max(im2(:)),20);
    m=m>0;
    try
        [voxsize,center]=dcmDimCenter(dname);
        
        if size(im,4)>1
            d=im2;
            save([prefix,'_mean.mat'],'d','voxsize','center');
        end
        
    catch
        if size(im,4)>1 && ~strcmp(dname(end-3:end),'.mat')  % not a mat file
            save([prefix,'_mean.mat'],'im2');
        end
    end
    
    
    figure;imshow4(double(m),[],[1,size(m,3)]);
    
else
    m=roi_name;
end


a=double(ri_d1(dname,1));

if exist('inv_off','var') && inv_off>0
    a=-a;   % inversion_off in sequence; inversion manually here
end


a_max=max(a(:));
a_min=min(a(:));

a=(a-(a_min+a_max)/2)/(a_max-a_min)*360;

a_samebg=a*0;

a_indvbg=a*0;

for sl=1:size(a,3)
    pos=ind2subb(size(m(:,:,sl)),find(m(:,:,sl)>0));
    x=[pos,ones(size(pos,1),1),pos.^2];
    
    ix=(x'*x)\x';
    a_tmp=mean(a(:,:,sl,:),4);
    y=a_tmp(m(:,:,sl)>0);
    b_mean=ix*y;
    
    for k=1:size(a,4)
        a_tmp=a(:,:,sl,k);
        y=a_tmp(m(:,:,sl)>0);
        b_indvbg=ix*y;
        
        for i=1:size(a,1)
            for j=1:size(a,2)
                a_samebg(i,j,sl,k)=a(i,j,sl,k)-[i,j,1,i^2,j^2]*b_mean;
                if size(a,4)>1
                    a_indvbg(i,j,sl,k)=a(i,j,sl,k)-[i,j,1,i^2,j^2]*b_indvbg;
                end
            end
        end
    end
    
end
%%

    prefix=strtok2(prefix,'.');

if inv_off
    prefix=[prefix,'_inv'];
end

try
    [voxsize,center]=dcmDimCenter(dname);
    
    if size(a,4)>1
        d=a_indvbg;
        save([prefix,'_detrend_indvbg.mat'],'d','voxsize','center');
        
        d=a_samebg;
        save([prefix,'_detrend_samebg.mat'],'d','voxsize','center');
        
        d=mean(a_samebg,4);
        save([prefix,'_detrend_mean.mat'],'d','voxsize','center');
    else
        d=a_samebg;
        save([prefix,'_detrend.mat'],'d','voxsize','center');
    end
    
catch
    if size(a,4)>1
        save([prefix,'_detrend_indvbg.mat'],'a_indvbg');
        
        save([prefix,'_detrend_samebg.mat'],'a_samebg');
        
        ma=mean(a_samebg,4);
        save([prefix,'_detrend_mean.mat'],'ma');
    else
        save([prefix,'_detrend.mat'],'a_samebg');
    end
    
end




