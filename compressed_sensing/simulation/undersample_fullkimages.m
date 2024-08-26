[a,info]=BrikLoad('synthesize_kdata_nois0.040_bold0.020+orig');

WriteBrikEZ(a(:,:,:,1:4:end-3),info,'undersampled from synthesize_kdata_nois0.040_bold0.020+orig','synthesize_kdata_nois0.040_bold0.020_R4','');

b=load('refg_0_10_30_6cy_TR2.0');
c=b(1:4:end-3);
fid=fopen('refg_0_10_30_6cy_TR2.0_R4','w');
for i=1:length(c)
    
    fprintf(fid,'%f\n',c(i));
end
fclose(fid);
    

correlateAnalysis('synthesize_kdata_nois0.040_bold0.020_R4+orig','refg_0_10_30_6cy_TR2.0_R4','cc_synthesize_kdata_nois0.040_bold0.020_R4');
