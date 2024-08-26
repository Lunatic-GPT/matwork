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

for i=1:length(prefix)
     for k=1:5
       %{  
         if k==3 && i==2
             fname = 'modCS_gauss_nois0.040000_xfmW0.100';
         else
         %}
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
         
          %  fname = sprintf('synthesize_kdata_nois%4.3f_bold%4.3f',nois(k),bold_scl);
       %     if exist(['cc_',fname,'+orig.HEAD'],'file')    
             cc=BrikLoadf(['cc_',fname,'+orig.HEAD[1]']);
             
             for l=1:length(thr)
                                  
              tmp=cc>thr(l) &mask>0;
              tmp2=(cc)>thr(l)&mask==0;
              tmp4=tmp2&mask2;
              %{
              [tmp3,clus]=clusterize2(cc>=thr(l),10);
              m=zeros(size(cc));
              
              if ~isempty(clus)
              for i=1:size(clus{1},1)
                  m(clus{1}(i,:))=1;
              end
              end
              if length(clus)>1
              for i=1:size(clus{2},1)
                  m(clus{2}(i,:))=1;
              end
              end
              
              fdetf2(k,l)=length(find(tmp3>0&mask==0))/nact;
              %}
              a=length(find(tmp));
              detf(i,k,l)=a/nact;
              fdetf(i,k,l)=length(find(tmp2))/ninact;
              fdetf3(i,k,l)=length(find(tmp4>0))/length(find(mask==0&mask2>0));
              
             end
             
             tmp=(detf(i,k,1:end-1)+detf(i,k,2:end))/2.*(fdetf(i,k,2:end)-fdetf(i,k,1:end-1));
             areaf(i,k)=-sum(tmp(:));
             tmp=(detf(i,k,1:end-1)+detf(i,k,2:end))/2.*(fdetf3(i,k,2:end)-fdetf3(i,k,1:end-1));
             areaf2(i,k)=-sum(tmp(:));
        %    end                
     end
end
%%

 
sym={'>-r','s-b','o-g','x-k','x-r'};
tlt={'kt','mod-cs','spatial'};
for k=1:3
    figure;
    for l=1:5
     hold on;plot(squeeze(fdetf(k,l,:)),squeeze(detf(k,l,:)),sym{l},'MarkerSize',8); 
    end
    lgstr=pattern;
    lgstr{end+1}='full';
    legend(lgstr);
    title(tlt{k});
    
end




































