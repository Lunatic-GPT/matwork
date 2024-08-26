function T1T2map_RAREVTR_Hong_dcm()
% T1map_sdt(fid_prefix,[t1t2_0])
% fid_prefix.
% t1t2_0: initial values for T1 and T2. default: 2, 0.045
% fit_mean: fit averaged data; default true;
% only tested for single slice 
if ~exist('stretch_exp','var') || isempty(t1t2)
    t1t2_0 = [2,0.045];
end

if ~exist('fit_mean','var')
    fit_mean=true;
end
dir_str=dir('*.dcm');

tr=[];
img=[];
pos=[];
for i=1:length(dir_str)
    h=dicominfo(dir_str(i).name);
    tmp=dicomread(dir_str(i).name);
    img=cat(3,img,tmp);
    tr(i)=h.RepetitionTime;
    pos(i)=h.SliceLocation;
end

tr2=unique(tr);
pos2=unique(pos);

tr2=sort(tr2,'ascend');

pos2=sort(pos2,'ascend');

if size(img,3)~=length(tr2)*length(pos2)
    error('Number of images wrong');
end

img2=zeros(size(img,1),size(img,2),length(pos2),length(tr2));

for i=1:size(img,3)
   img2(:,:,pos(i)==pos2,tr(i)==tr2)=img(:,:,i);
    
end

vd=tr2/1000;


sz=size(img2);

R1=zeros(sz(1:3));
res1=zeros(sz(1:3));


options=optimset('MaxIter',20,'Display','off');

tmp=img2(:,:,:,end);
m=tmp>0.02*max(tmp(:));

for i=1:sz(1)
    fprintf('%d/%d\n',i,sz(1));
    for j=1:sz(2)
        for k=1:sz(3)
          if m(i,j,k)==0
              continue;
          end
            
            
            y = squeeze(img2(i,j,k,:));
        
            [beta,r]=lsqcurvefit(@exp_decay2,[-max(y),t1t2_0(1),max(y)],vd(:),y(:),[],[],options);
   %        disp(toc); 
            if ~any(isnan(beta)) 
             ss = sum((y-mean(y)).^2);   
             R1(i,j,k) = 1/beta(2);
             res1(i,j,k) = 1-sum(r.^2)/ss;
            end
    
        end
           
    end
end

  save('R1','img2','R1');


function y=exp_decay2(b,x)
% y=b1*exp(-x/b2)+b3;
% %4.3f*exp(-x/%4.3f)+%4.3f
y=b(1)*exp(-x/b(2))+b(3);


