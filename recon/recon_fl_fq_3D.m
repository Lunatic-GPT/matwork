function recon_fl_fq_3D(fmb)

%fmb='meas_MID42_fl_fqZ_AF3_FID10324.dat';
%fsb='meas_MID44_fl_fqZ_SB_FID10326.dat';

% GRAPPA recon
% 
% fmb='meas_MID72_fl_fq_mb_VENC2_640_TEopt_BWopt_FID12373.dat';
% fsb='meas_MID70_fl_fq_sb_VENC2_640_TEopt_BWopt_FID12371.dat';
% 
% fmb='meas_MID85_fl_fq_mb_VENC2_640_TEopt_BWopt_FID12486.dat';
% fsb='meas_MID90_fl_fq_sb_VENC2_640_TEopt_BWopt_FID12491.dat';
% 
% fmb='meas_MID59_fl_fq_mb_VENC4_640_TEopt_BWopt_FID12925.dat';
% fsb='meas_MID61_fl_fq_sbRef_VENC4_640_TEopt_BWopt_FID12927.dat';

%{
fmb='meas_MID20_fl_fq_retroZ_MB_sep2_FID12512.dat';
fsb='meas_MID22_fl_fq_retroZ_SB_sep2_FID12514.dat';

fmb='meas_MID21_fl_fq_retroZ_MB_sep2_smallFOV_FID12513.dat';
fsb='meas_MID23_fl_fq_retroZ_SB_sep2_smallFOV_FID12515.dat';

fmb='meas_MID33_fl_fq_mb_VENC2_640_TEopt_BWopt_FID12525.dat';
fsb='meas_MID36_fl_fq_sb_VENC2_640_TEopt_BWopt_FID12528.dat';
%}
%fmb='meas_MID24_fl_fq_retroZ_mb_SINC_FID13783.dat';


[d,lin,par]=readMeasDat(fmb,inf,0,true);
%fmb='meas_MID78_fl_fqZ_MB_TE14_5ms_FA35_FID11035.mat';
%fsb='meas_MID76_fl_fqZ_SBRef_TE13_2ms_FA35_FID11033.mat';

lin2=lin(1:64:end);
par2=par(1:64:end);
nlin=max(lin)+1;
npar=max(par)+1;
nro=size(d,1);
d=reshape(d,[nro,32,2,length(lin2)]);

d2=zeros(nro,32,2,nlin,npar);
for i=1:length(lin2)
    
    d2(:,:,:,lin2(i)+1,par2(i)+1)=d(:,:,:,i);
end

fd2=ifft1c(ifft1c(d2,4),5);


afd2=squeeze(sos(fd2,2));

afd2=permute(afd2,[1,3,4,2]);
prefix=strtok(fmb,'.');
save(prefix,'afd2');
