
tvWeight=[0.001,0.01];
ftWeight=[0,0.001,0.005,0.01];
nois=[0.01,0.04,0.1];

%% fully sampled; correlation analysis
bold_scl=0.02;
tv=[0,0.002];
ftw=[0,0.002];
pattern={'uniform','uniform_nc6','gauss'};

for i=1:2
    for j=1:2
      for k=1:length(pattern)
          if i==1 && j==1
              continue;
          end
          
       fname = sprintf('run_cs_ft_nois0.010_tv%4.3f_ft%4.3f%s',tv(i),ftw(j),pattern{k});
       if exist([fname,'+orig.HEAD'],'file') && ~exist(['cc_',fname,'+orig.HEAD'],'file')
        correlateAnalysis([fname,'+orig'],'refg_0_10_30_6cy_TR2.0',['cc_',fname]);
       end
      end
    end
end
%% ROC curves

thr=[linspace(0,0.3,20),0.4,0.6,0.8,1];
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
detf=zeros(3,4,length(thr)); %3 regularization, 3 patterns
fdetf=zeros(3,4,length(thr)); 

fdetf2=zeros(3,4,length(thr));  %clusterized and normalized by the total number of activated voxels.


fdetf3=zeros(3,4,length(thr));% same as fdetf but only consider voxles within mask2
areaf=zeros(3,4);
areaf2=zeros(3,4);

tv_ftw=[0,0.002;0.002,0;0.002,0.002];

for i=1
     for k=1:5
         
         if k<4
          fname=sprintf('run_cs_ft_nois0.010_tv%4.3f_ft%4.3f%s',tv_ftw(i,1),tv_ftw(i,2),pattern{k});
         else
          fname='synthesize_kdata_nois0.010_bold0.020';
         end
         
          %  fname = sprintf('synthesize_kdata_nois%4.3f_bold%4.3f',nois(k),bold_scl);
            if exist(['cc_',fname,'+orig.HEAD'],'file')    
             cc=BrikLoadf(['cc_',fname,'+orig.HEAD[1]']);
             
             for l=1:length(thr)
                                  
              tmp=cc>thr(l) &mask>0;
              tmp2=abs(cc)>thr(l)&mask==0;
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
            end                
     end
end
%%

 
sym={'>-r','s-b','o-g','x-k'};
tlt={'FT only','TV only','FT+TV'};
for k=1:3
    figure;
    for l=1:4   
     hold on;plot(squeeze(fdetf(k,l,:)),squeeze(detf(k,l,:)),sym{l},'MarkerSize',8); 
    end
    lgstr=pattern;
    lgstr{end+1}='full';
    legend(lgstr);
    title(tlt{k});
end


%% under sampled; correlation analysis

         
for i=1:length(tvWeight)
    for j=1:length(ftWeight)
        for k=1:length(nois) 
            fname = sprintf('run_cs_ft_nois%4.3f_tv%4.3f_ft%4.3f+orig',nois(k),tvWeight(i),ftWeight(j));
            if exist([fname,'.HEAD'],'file')  && ~exist(['cc_',fname,'.HEAD'],'file')
                
             correlateAnalysis(fname,'refg_0_10_30_6cy_TR2.0',['cc_',fname]);
            end
        end
    end
end

%% ROC curves

mask=zeros(64,64,9);
mask(19:25,45:50,4:6)=1;
mask(45:51,43:48,4:6)=1;

nact=length(find(mask>0));
ninact=length(mask(:))-nact;
thr=[linspace(0,0.3,10),0.4,0.6,0.8];

det=zeros(length(tvWeight),length(ftWeight),length(nois),length(thr));
fdet=zeros(length(tvWeight),length(ftWeight),length(nois),length(thr));
fdet2=zeros(length(tvWeight),length(ftWeight),length(nois),length(thr));
fdet3=zeros(length(tvWeight),length(ftWeight),length(nois),length(thr));

area1=zeros(length(tvWeight),length(ftWeight),length(nois));
area2=zeros(length(tvWeight),length(ftWeight),length(nois));

iplot=1;
for i=1:length(tvWeight)
    for j=1:length(ftWeight)
        for k=1:length(nois)
            fname = sprintf('run_cs_ft_nois%4.3f_tv%4.3f_ft%4.3f+orig',nois(k),tvWeight(i),ftWeight(j));
            if exist(['cc_',fname,'.HEAD'],'file')    
             cc=BrikLoadf(['cc_',fname,'.HEAD[1]']);
             
             for l=1:length(thr)
              tmp=cc>thr(l) &mask>0;
              tmp2=abs(cc)>thr(l)&mask==0;
              tmp4=tmp2&mask2>0;
         %     [tmp3,clus]=clusterize2(cc>=thr(l),10);

              a=length(find(tmp));
              det(i,j,k,l)=a/nact;
              fdet(i,j,k,l)=length(find(tmp2))/ninact;
              
              fdet3(i,j,k,l)=length(find(tmp4>0))/length(find(mask==0&mask2>0));
        %      fdet2(i,j,k,l)=length(find(tmp3>0&mask==0))/nact;
             end
             tmp=(det(i,j,k,1:end-1)+det(i,j,k,2:end))/2.*diff(fdet(i,j,k,:));
             area1(i,j,k)=-sum(tmp(:));
             
             tmp=(det(i,j,k,1:end-1)+det(i,j,k,2:end))/2.*diff(fdet3(i,j,k,:));
             area2(i,j,k)=-sum(tmp(:));
            end                     
        end
    end
end
%% ROC curves
ntv=length(tvWeight);

figure;
nft=length(ftWeight);
sym={'>-r','s-b','o-g','x-k'};
  for i=1:ntv
    for j=1:nft
          subplot(ntv,nft,(i-1)*nft+j);
          for k=1:length(nois)
        
           x=squeeze(fdet(i,j,k,:));
           y=squeeze(det(i,j,k,:));
           
           hold on;plot(x,y,sym{k});   
           title(sprintf('tvW=%4.3f; ftW=%4.3f',tvWeight(i),ftWeight(j)));
          end
          
          xlim([0,0.5]);
    end
  end
         

  
%% Area under ROC curves
figure;
ntv=length(tvWeight);
nft=length(ftWeight);
sym={'>-r','s-b','o-g','x-k'};
  for i=1:ntv
    for j=1:nft
          subplot(ntv,nft,(i-1)*nft+j);
          for k=1:length(nois)
        
           y=[areaf2,squeeze(area2(i,j,:))];
           
           h =bar(y);
           set(h(1),'FaceColor',[0.1,0.1,0.1])
           set(h(2),'FaceColor',[0.8,0.8,0.8])
           title(sprintf('TV=%4.3f; FT=%4.3f',tvWeight(i),ftWeight(j)));
          end      
          xlim([0.5,3.5]);
    end
  end

































