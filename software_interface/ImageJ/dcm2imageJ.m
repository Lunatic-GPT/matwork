d=read_afni_sdt_images('22');

d2=floor((d)/600*255);
d2(d2>254)=254;


roi=load('PVS_final_WithHmethod.mat');

d2(roi.roi>0)=255;
writeanalyze(d2,'22_Hmethod',[1,1,1],'uint8');


%%

