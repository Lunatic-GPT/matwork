nois_array=[0.04,0.10];

for j=1
    
 nois=nois_array(j);
 tmpname = sprintf('synthesize_kdata_mask2_nois%4.3f_bold%4.3f.mat',nois,0.02);
 tmp=load(tmpname);
 z2=tmp.z3;


 pat={'uniform','uniform_nc6','gauss'};
 for i=1:3

   m=get_mask(pat{i});
        
   
   run_cs_ft_mask2(nois,0.000,0.001,m,pat{i});
   
          
end

end