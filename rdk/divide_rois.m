retino_wdgrng('wedge_vr_ave_rs_remote_blur+orig','ring_vr_ave_rs_remote_blur+orig','retino_remote',16)
!3dcalc -prefix mask_ecc -a retino_remote+orig'[2]' -b retino_remote+orig'[3]' -expr '(1+step(a-180))*ispositive(b-0.24)' 
maskcalc(mr{7},ml{7},'or(a,b)',1,'mask_MTLoc');
separateMask('mask_MTLoc+orig','ap','mask_MTLoc_ap');

separateMask('V5m+tlrc',3,'ap','V5m_ap');
separateMask('V5u+tlrc',3,'ap','V5u_ap');

maskcalc('mask_rdk_REML+orig','equals(a,1)+equals(a,2)',1,'mask_MTLoc');