function dcm2afni_siemens(dname,mtrx)
% dcm2afni_siemens(dname,mtrx)
tic;

     str = dir(fullfile(dname,'*.IMA'));

     a=dicominfo(fullfile(dname,str(1).name));
%mtrx=[a.Width,a.Height];

     nt = length(str);
     
     if nt==0
         str=dir(fullfile(dname,'*.dcm'));
         nt = length(str);
     
     end
     
     brikdata=zeros([mtrx,nt]);
    for i=1:nt
        
        img = dicomread(fullfile(dname,str(i).name));
        tmp=zeros(mtrx);
        nx=size(img,2)/mtrx(2);
        ny=size(img,1)/mtrx(1);
        
        for j=1:mtrx(3)
            x=mod(j-1,nx);
            y=floor((j-1)/nx);
           tmp(:,:,j)= img(y*mtrx(1)+1:(y+1)*mtrx(1),x*mtrx(2)+1:(x+1)*mtrx(2));
        end
        brikdata(:,:,:,i) = tmp;
    end
    
        


       
       write_afni(brikdata,dname);
         
       
     
%% save the average time series for each trial type.    

    disp([mfilename ' finish in ', num2str(toc), ' s']);