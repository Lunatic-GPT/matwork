pattern={'uniform','uniform_nc6','gauss'};

nois=0.04;
%prefix={'run_cs_ft_nois0.040_tv0.000_ft0.004%s','modCS_%s_nois0.040000_xfmW0.050','run_cs_wavelet_nois0.040_tv0.001_ft0.001%s'};
 prefix = {'ktFOCUSS_nih_nois0.020_%s','modCS_nih_%s_nois0.020000_bold0.015_xfmW0.050','run_cs_wavelet_nois0.040_tv0.001_ft0.001%s'};   
for i=1:length(prefix)
   for j=1:length(pattern)
          
       fname = sprintf(prefix{i},pattern{j});
       if exist([fname,'+orig.HEAD'],'file') && ~exist(['cc_',fname,'+orig.HEAD'],'file')
        correlateAnalysis([fname,'+orig'],'refg_0_10_30_6cy_TR2.0',['cc_',fname]);
       end
   end
end



%% ROC curves

thr=[linspace(0,0.1,100),linspace(0.101,0.3,30),0.4,0.6,0.8,1];
nt=240;
mask=zeros(64,64,9);
mask2=zeros(64,64,9);
mask(19:25,45:50,4:6)=1;
mask(45:51,43:48,4:6)=1;

mask2(17:27,43:52,3:7)=1;
mask2(43:53,41:50,3:7)=1;

nact=length(find(mask>0));
ninact=length(mask(:))-nact;
%thr=linspace(0,1,20);
detf=zeros(3,5,length(thr)); %3 regularization, 3 patterns
fdetf=zeros(3,5,length(thr)); 

fdetf2=zeros(3,5,length(thr));  %clusterized and normalized by the total number of activated voxels.

fdetf3=zeros(3,5,length(thr));% same as fdetf but only consider voxles within mask2
areaf=zeros(3,5);
areaf2=zeros(3,5);

tv_ftw=[0,0.002;0.002,0;0.002,0.002];
ninact=length(find(mask==0));

pval=zeros(ninact,length(prefix),5);

for i=1:length(prefix)
     for k=1:5

         if k==2 && i==2
             fname = 'modCS_nih_uniform_nc6_nois0.020000_bold0.015_xfmW0.100';
             
         elseif k==1 && i==1
             fname = 'run_cs_ft_nois0.040_tv0.000_ft0.002uniform';
         elseif k<4
          fname=sprintf(prefix{i},pattern{k});
          
         elseif k==4
          fname='synthesize_kdata_nih_nois0.020_bold0.015';
         else
             fname = 'synthesize_kdata_nih_nois0.020_bold0.015_R4';
         end
         
  
             t=BrikLoadf(['cc_',fname,'+orig.HEAD[2]']);
            ttmp=t(mask==0);
            disp([i,k]);
             for j=1:ninact
               pval(j,i,k)=tTest(238,ttmp(j));
             end
        
     end
end
%%












