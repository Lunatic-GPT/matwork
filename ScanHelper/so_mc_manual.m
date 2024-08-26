function so_mc_manual(pvsscan,navscan)
% pvsscan: the scan # for pvs image;
% navscan: the scan # for navigators;
% To be run in ubuntu; 
% The pvsscan should already be converted to nii file
% by dcm2pvs

if isa(navscan,'char')
  navscan=str2double(navscan);
end

if isa(pvsscan,'char')
  pvsscan=str2double(pvsscan);
end


fprintf('\nUsing pvsscan %d\n',pvsscan)
ref=scan_dir_name(pvsscan);




%fpos_base={'SlicePosition//SlicePosition1_base.txt','SlicePosition//SlicePosition2_base.txt','SlicePosition//SlicePosition3_base.txt'};
%fpos={'SlicePosition//SlicePosition1.txt','SlicePosition//SlicePosition2.txt','SlicePosition//SlicePosition3.txt'};


fpos_base={'SlicePosition//SlicePosition3_base.txt'};
fpos={'SlicePosition//SlicePosition3.txt'};

dname=scan_dir_name(navscan);

%if ~exist([dname,'.nii.gz'],'file')  %has not been processed
    tic;

  ScanSliceAdjust(ref,dname,fpos_base,fpos);
   
    
   
    toc;
%end



function plot_motion_afni(d)


lbl={'IS','RL','AP','IS','RL','AP'};

for i=1:3
    subplot(3,2,i*2-1);
    hold off;
    plot(d(:,3+i),'o-');
    hold on;
    plot([0.5,size(d,1)+0.5],[0,0],'--');
    ylabel([lbl{3+i},' (mm)']);
    
    if i==1
        title('Tranlation');
    end
    xlim([0.5,size(d,1)+0.5]);
end

for i=1:3
    subplot(3,2,i*2);
    hold off;
    plot(d(:,i),'o-');
    hold on;
    plot([0.5,size(d,1)+0.5],[0,0],'--');
    ylabel([lbl{i},' (deg)']);
    if i==1
        title('Rotation');
    end
    xlim([0.5,0.5+size(d,1)]);
    
end

function res=scan_dir_name(i)

a=dir(sprintf('*_%04d',i));


if strcmp(a(1).name(1:4),'mask') %make sure it is not a PVS mask
    res=a(2).name;
else
    res=a(1).name;
end




