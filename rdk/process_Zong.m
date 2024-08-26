cd /media/SS_FMRI1/MRI_Analysis/Hua/rdk_12_04_08/afni

%%
cmd = '3dAutomask -prefix brainmask rdk1_vr+orig';
unix(cmd);

rois = {'V1','V2d','V2v','V3d','V3v','V4','MTLoc'};
for i=1:length(rois)
    rois_l{i} = ['../../rois/',rois{i},'.lh.1D.roi'];
    rois_r{i} = ['../../rois/',rois{i},'.rh.1D.roi'];
end
sroi2v('brainmask+orig',rois_l,'Hua','lh',1);
sroi2v('brainmask+orig',rois_r,'Hua','rh',1);

%%
tps = [6:20,35:40,55:60,75:80,95:100,115:120,135:140,155:160];
tps = tps - 5;
for i=1:10
  fnames{i} = sprintf('rdk%d_dtr_norm+orig',i);
  cmd = sprintf('3dFourier -prefix rdk%d_dtr_norm_lp0p1 -lowpass 0.1 %s',i,fnames{i});
  unix(cmd);
end
tsstd(fnames,tps,'rdkAll_std_baseline');
tssort_bin(fnames,'stim/stim_cohall.1D',20,5,0,'tssorted');


for i =1:10
   fnames{i} = sprintf('rdk%d_dtr_norm_lp0p1+orig',i);
end
tsave_rdk(fnames,'stim/stim_cohall.1D',20,5,'tsave_lp0p1');
tssort_bin(fnames,'stim/stim_cohall.1D',20,5,0,'tssorted_lp0p1');

ref_gamma_conv(0,5,15,20,'ref_gamma_20trials.1D')
for i=1:4
correlateAnalysis(['tssorted_lp0p1_coh',num2str(i),'+orig'],'ref_gamma_20trials.1D',['cc_lp0p1_coh',num2str(i)]);
correlateAnalysis(['tssorted_coh',num2str(i),'+orig'],'ref_gamma_20trials.1D',['cc_coh',num2str(i)]);
end

%% 
load ../rawData/12-04_15-44_RDK_fMRI_Hua.mat
ncoh = length(params.coh);
coh = zero(1,ncoh+1); 
coh(2:ncoh+1) = params.coh;

[lroi,nl]=mask_str('brainmask_sroi2v_lh+orig');
[rroi,nr]=mask_str('brainmask_sroi2v_rh+orig');
nrois = length(lroi);
thr_list = [0.15,0.2,0.3];

for i=1:ncoh
    ts{i} = sprintf('tsave_lp0p1_coh%d+orig',i);
    ts_sort{i} = sprintf('tssorted_coh%d+orig',i);
    cc{i} = sprintf('cc_lp0p1_coh%d+orig[0]',i);
    cc{i} = sprintf('cc_coh%d+orig[0]',i);
end

symb = {'sr','*g','>b','+k'};
figure;
for it=1:length(thr_list)
    cc_all = [];
   for i=1:ncoh
    cc_all = sprintf('%sstep(%c-%3.2f)*',cc_all,'a'+i-1,thr_list(it));
   end 
    cc1 = sprintf('step(a-%3.2f)*',thr_list(it));

   for i=1:nrois
   subplot(3,ceil(nrois/3),i);   
   nv = zeros(1,ncoh+1);
   for j=1:ncoh+1
     if j==1
        expr = [cc_all,sprintf('(%c+%c)','a'+ncoh,'a'+ncoh+1)];
        nv(1) = nv_roi('brainmask+orig',maskcalc(cc{1:ncoh},rroi{i},lroi{i},expr)); 
     else
        expr = [cc1,'(b+c)'];
        nv(j)=nv_roi('brainmask+orig',maskcalc(cc{j-1},rroi{i},lroi{i},expr));
     end
   end
   hold on;
   plot(coh,nv,symb{it});
   end
end

for i =1:nrois
    subplot(3,ceil(nrois/3),i);
  plot_labels([nl{i},'+',nr{i}],'coherence','number of voxels',[],16);
end
