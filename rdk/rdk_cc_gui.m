function rdk_cc_gui
 sid = 2;
 if sid == 1 || sid == 3
   ref_gamma_conv(0,5,15,20,'ref_gamma_20.1D');
   for i=1:4
    flist{i} = sprintf('ts_sort_coh%d+orig',i);
   end
   correlateAnalysis(flist,'ref_gamma_20.1D','cc_gamma');
 elseif sid == 2
   ref_gamma_conv(0,5,15,20,'ref_gamma_20.1D');  % Peng
   for i=1:4
    flist = sprintf('ts_sort_coh%d+orig',i);
    correlateAnalysis(flist,'ref_gamma_20.1D',sprintf('cc_gamma_coh%d',i));
   end
   !3dTcat -prefix cc_gamma cc_gamma_coh1+orig cc_gamma_coh2+orig cc_gamma_coh3+orig cc_gamma_coh4+orig
 else
   ref_gamma_conv(0,5,15,32,'ref_gamma_32.1D');  % K004
   for i=1:3
    flist = sprintf('ts_sort_coh%d+orig',i);
    correlateAnalysis(flist,'ref_gamma_32.1D',sprintf('cc_gamma_coh%d',i));
   end
   !3dTcat -prefix cc_gamma cc_gamma_coh1+orig cc_gamma_coh2+orig cc_gamma_coh3+orig
end
