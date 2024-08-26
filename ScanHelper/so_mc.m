function so_mc(pvsscan,navscan)
% which pvsscan to use; first one found if 0;
% which navscan to use; latest one if 0;
%
if ~exist('pvsscan','var')
    pvsscan=0;
end

if ~exist('navscan','var')
    navscan=0;
end

tic;
fprintf('Waiting for pvs scan 0000');

% if not specified, then search
if pvsscan==0
  while 1
    if scan_exist(pvsscan) && ispvs(pvsscan)
        break;
    end
    
    if scan_exist(pvsscan+1)
        pvsscan=pvsscan+1;
    end
    fprintf('\b\b\b\b%04d',pvsscan);
    pause(1);
  end
end

fprintf('\nUsing pvsscan %d\n',pvsscan);

ref=scan_dir_name(pvsscan);

if ~exist([ref,'.nii.gz'],'file')
    dcm2nii(ref);
end


motion_par=[];


figure(11);

fprintf('Waiting for nav scan %04d',navscan);
while 1
    
    fpos_base={'SlicePosition//SlicePosition1_base.txt','SlicePosition//SlicePosition2_base.txt','SlicePosition//SlicePosition3_base.txt'};
    fpos={'SlicePosition//SlicePosition1.txt','SlicePosition//SlicePosition2.txt','SlicePosition//SlicePosition3.txt'};
    
    if scan_exist(navscan) && isnav(navscan)
     
        dname=scan_dir_name(navscan);
        
        if ~exist([dname,'.nii.gz'],'file')  %has not been processed
            tic;
            if isempty(motion_par)
                motion_par= ScanSliceAdjust(ref,dname,fpos_base,fpos);
            else
                motion_par(end+1,:)= ScanSliceAdjust(ref,dname,fpos_base,fpos);
            end
            
            plot_motion_afni(motion_par);
            toc;
        end
        
    end
    pause(1);
    
    if scan_exist(navscan+1)
        navscan=navscan+1;
    end
    fprintf('\b\b\b\b%04d',navscan);
end

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

function res=scan_exist(i)

a=dir(sprintf('*_%04d',i));

res= ~isempty(a);

function res=scan_dir_name(i)

a=dir(sprintf('*_%04d',i));


if strcmp(a(1).name(1:4),'mask') %make sure it is not a PVS mask
    res=a(2).name;
else
    res=a(1).name;
end


%% assume the file name convention will not change after files removed from the output folder

function resLast=last_scan_id(nav_pvs)
% the pvsscan'th pvs scan is ready.

resLast=0;
for i=15:-1:1
    disp(i);
    if ~scan_exist(i)
        continue;
    end
    
    if nav_pvs==0
        tmp=isnav(i);
    else
        tmp=ispvs(i);
    end
    
    if tmp
        resLast=i;
        break;
    end
    
    
    
end


function res=ispvs(i)
%% check whether the scan is a pvs scan and has all the files ready
res=false;

dname=scan_dir_name(i);

a=readdPar(dname,'ProtocolName');
if length(a)>=17 && strcmp(a(1:17),'tse_vfl_pss_FNRecon_matchR21')
    
    dir_str=dir(fullfile(dname,'*.dcm'));
    
    if length(dir_str)==248
        res=true;
    end
    
end





function res=isnav(i)
%% check whether the scan is a pvs scan and has all the files ready
res=false;

dname=scan_dir_name(i);

a=readdPar(dname,'ProtocolName');
if length(a)>=15 && strcmp(a(1:15),'tse_vfl_pss_Nav')
    
    dir_str=dir(fullfile(dname,'*.dcm'));
    if isempty(dir_str)
     dir_str=dir(fullfile(dname,'*.IMA'));   
    end
    
    if length(dir_str)==52
        res=true;
    end
    
end


function res3=readdPar(dname,par,allfile)
%res=readdPar(dname,par,allfile)

d=dir(dname);
d=d(3:end);

warning off;
if ~exist('allfile','var')  || ~allfile
    
    for i=1:length(d)
        if isdir(fullfile(dname,d(i).name))
            continue;
        end
        try
            in=dicominfo(fullfile(dname,d(i).name));
            
            break;
        catch
            continue;
        end
    end
    if exist('par','var')
        res3=getfield(in,par);
    else
        
        disp(in);
    end
    
else
    
    res3={};
    for i=1:length(d)
        %   disp(i);
        in=dicominfo(fullfile(dname,d(i).name));
        try
            res3{i}=getfield(in,par);
        catch
            res3{i}='';
        end
        
    end
    
    
    if isa(res3{1},'char')
        return;
    end
    try
        res3=cell2mat(res3);
    catch
        
    end
    
end

