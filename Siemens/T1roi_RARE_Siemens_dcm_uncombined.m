function T1=T1roi_RARE_Siemens_dcm_uncombined(dname,mask,channels,mask_label)
% T1roi_RARE_Siemens_dcm_uncombined(dname,mask,channels,mask_label)
% use uncombined data
% dname is an array (doubl or cell) of dicom directories. 
% channels: one based.  Assume 32 channels
% mask can be matrix, striing cell array, or string


%dname={'9','10','11','12','13','14'};

if isa(dname,'double')
   
    for i=1:length(dname)
        dname2{i}=num2str(dname(i));
    end
    dname=dname2;
end
    
a=[];
for i=1:length(dname)

    tmp=ri(dname{i},32);
    
    a(:,:,:,i)=tmp(:,:,channels);
    vd(i)=readdPar(dname{i},'RepetitionTime');
end

a=rms(a,3);

[vd,ind]=sort(vd);
a=a(:,:,1,ind);

if isa(mask,'str') || isa(mask,'cell')
    
  mask=str2cell(mask);
else
  mask={mask};
end

for i=1:length(mask)
if isa(mask{i},'char')
    m=load(mask{i});
    m=m.roi;
else
    m=mask{i};
end


sz=size(a);

R1=zeros(sz(1:3));
res1=zeros(sz(1:3));

options=optimset('MaxIter',20,'Display','off');


            
            y=mean_roi(a,m);
            [beta,r]=lsqcurvefit(@exp_decay2,[-max(y),mean(vd),max(y)],vd(:),y(:),[],[],options);
   %        disp(toc); 
             R1 = 1/beta(2);
            figure;plot(vd,y,'o');
           hold on;plot(vd,exp_decay2(beta,vd),'r-');
           
           
if isa(mask{i},'char')
           fprintf('%s: T1 = %f\n',mask{i},1/R1);
           title( sprintf('%s: T1 = %f\n',mask{i},1/R1));
else
       fprintf('%s: T1 = %f\n',mask_label,1/R1);
       title(sprintf('%s: T1 = %f\n',mask_label,1/R1));
end
       T1(i)=1/R1;     
end
          
           
             
function y=exp_decay2(b,x)
% y=b1*exp(-x/b2)+b3;
% %4.3f*exp(-x/%4.3f)+%4.3f
y=b(1)*exp(-x/b(2))+b(3);


