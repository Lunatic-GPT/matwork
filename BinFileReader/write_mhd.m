%roi=load('../PVSLength_testzb_Curv_Dv_Comb_Dv.mat');


roi=permute(roi.c,[2,1,3]);

f=fopen('PVS_final_WithHMethod_corrected2_xz.raw','wb');
%f=fopen(fname,'rb');

d=fwrite(f,roi,'integer*2');

fclose(f);

%copyfile PVS_final_WithHmethod.mhd PVS_final_WithHMethod_corrected2_xz.mdh

