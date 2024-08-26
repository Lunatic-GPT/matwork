%function recon_fl_fq_mb_grappa(fmb,fsb,fnois,interp_factor,do_retro,nphase)
% if fsb is an cell array; the position should be in ascending order, i.e. from foot
% to head.

remove_os=true;

dnois=readMeasDat(fnois,32,0,remove_os);

fsb=str2cell(fsb);


prefix_mb=strtok(fmb,'.');
prefix_sb=strtok(fsb{1},'.');

prefix_mb_dir=fullfile(prefix_mb,prefix_mb);

prefix_sb_dir=fullfile(prefix_sb,prefix_sb);

%fmb='meas_MID28_fl_fq_retroZ_mb_Tone_FID13919.dat';

if ~exist([prefix_mb_dir,'.mat'],'file')
    
    [dmb,lin_mb,par_mb,sl_mb,ushSet,timeStamp,freePara]=readMeasDat(fmb,inf,0,true);
    save([prefix_mb_dir,'.mat'],'dmb','lin_mb','par_mb','sl_mb','ushSet','freePara');
    
else
    load([prefix_mb_dir,'.mat'],'lin_mb','par_mb','sl_mb','ushSet','freePara');
    
    dmb=zeros(1,32,length(lin_mb)/32);
end

%[dmb,lin_mb,par_mb,sl_mb,ushSet,freePara]=readMeasDat(fmb,inf,0,true);
%fmb='meas_MID78_fl_fqZ_MB_TE14_5ms_FA35_FID11035.mat';
%fsb='meas_MID76_fl_fqZ_SBRef_TE13_2ms_FA35_FID11033.mat';
%%


if length(fsb)==1
    nsl=readsPar([prefix_sb_dir,'.pro'],'lConc');
else
    nsl=length(fsb);
end

af=readsPar([prefix_mb_dir,'.pro'],'alFree[6]');
nave=readsPar([prefix_mb_dir,'.pro'],'lAverages');
nave_sb=readsPar([prefix_sb_dir,'.pro'],'lAverages');

seg=readsPar([prefix_mb_dir,'.pro'],'lSegments');

%% sort physiology data

nvenc=max(ushSet)+1;

if ~do_retro
    nphase=1;
end

try
    phaseStab=readsPar([prefix_mb_dir,'.pro'],'ucPhaseStabilize');
    
    if strcmp(phaseStab{1},'0x1')
        phaseStab=true;
    else
        phaseStab=false;
    end
catch
    phaseStab=false;
end


if  phaseStab
    dmb2=reorder_fl_fq_data_phaseStabOn(dmb,lin_mb,prefix_mb_dir,nvenc,do_retro,freePara(:,4),nphase,false);  
else
    dmb2=reorder_fl_fq_data(dmb,lin_mb,prefix_mb_dir,nvenc,do_retro,freePara(:,4),nphase);
end
