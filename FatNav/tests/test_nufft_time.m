%{
% results
nufft: adjoin
nblock  imSize        time    memory      time/channel
1        half         4.5     0.246        0.025
2        half         4.6     0.123
1        full         5.3     0.2463       0.05
2        full         6.14    

nufft:
         ksize
1          full       5.94    0.245
1          half       2.93    0.123
2          full       5.8     0.123


%}
%%

ktmp=k2_center(1:end,:);

  nufft=NUFFT3D(ktmp.*repmat(imSize./imSize_lowres,[size(ktmp,1),1]),1,[0,0,0],round(imSize_lowres/4),1,1);
     tloc = tic;
  im = nufft'*Data(nskip*nro+ind_center(1:end),:);
  
  toc(tloc);
  
  
  %%
  
  
ktmp=k2_center(1:end,:);
   nufft=NUFFT3D(ktmp.*repmat(imSize./imSize_lowres,[size(ktmp,1),1]),1,[0,0,0],imSize_lowres,2,1);
     tloc = tic;
  tmp = nufft*im;
  
    toc(tloc);
