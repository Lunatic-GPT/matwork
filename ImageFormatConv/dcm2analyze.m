function pos=dcm2analyze(dname,sl_first,mosaic)
%dcm2analyze(dname,sl_first,mosaic)
% mosaic: [rows, cols] in the original image.

tic;

if ~exist('sl_first','var')
    sl_first=true;
end

tr=[];
img=[];
pos=[];

dir_str=dir2([dname,'/*']);
% if isempty(dir_str)
%  dir_str=dir([dname,'/*.IMA']);
% end    

zerobase=0;

for i=1:length(dir_str)
    
    
   if mod(i,10)==0
       disp(i);
   end
    
   if ~exist('ns','var')
    h=dicominfo(fullfile(dname,dir_str(i).name));
        
    if sl_first
        tmp=h.SliceLocation;
        if ~any(tmp==pos)
         pos(i)=tmp;
        else
            if nargout>0
                return;  %just get slice positions.
            end
            ns=length(pos);
        end
    else
         tmp=h.AcquisitionNumber;
    
       if ~any(tmp==tr)
            tr(i)=tmp;
       else
             ns=length(dir_str)/length(tr);
       end
   
   end
   end  
    img = dicomread(fullfile(dname,dir_str(i).name));
        
        if i==1
            dinfo = dicominfo(fullfile(dname,dir_str(i).name));
            brikdata = zeros([size(img),length(dir_str)]);
            vsize(1:2)=dinfo.PixelSpacing;
            vsize(3)=dinfo.SliceThickness;
            vsize(4)=1;
        end
        
        if ~isempty(brikdata)&& size(img,2)==size(brikdata,1) && size(img,1)==size(brikdata,2) && size(img,1)~=size(img,2)
              brikdata(:,:,i)=img';   %for some dicom images, the dimension switches.
            else
                brikdata(:,:,i)=img;
        end
            
        
end

if exist('mosaic','var')  && mosaic

    msize=readdPar(dname,'Private_0051_100b');   
    [msize1,msize2]=strtok(msize,'*');
    msize1=str2num(msize1);
    msize2=str2num(msize2(2:end));
    

brikdata=reshape(brikdata,[msize1,size(brikdata,1)/msize1,msize2,size(brikdata,2)/msize2,size(brikdata,3)]);
    
brikdata=permute(brikdata,[3,1,4,2,5]);

sz=size(brikdata);
brikdata=reshape(brikdata,[sz(1:2),prod(sz(3:4)),sz(5)]);
%save(dname,'brikdata');
     
%% save the average time series for each trial type.    

writeanalyze(brikdata,dname,vsize,'int16');
    disp([mfilename ' finish in ', num2str(toc), ' s']);
return;
end

if ~exist('ns','var') 
    if sl_first
        ns=size(brikdata,3);
    else
       ns=1;
    end
end

if sl_first
d=reshape(brikdata,[size(brikdata,1),size(brikdata,2),ns,size(brikdata,3)/ns]);
else 
 d=reshape(brikdata,[size(brikdata,1),size(brikdata,2),size(brikdata,3)/ns,ns]);
 d=permute(d,[1,2,4,3]);
end


writeanalyze(d,dname,vsize,'int16'); 

    disp([mfilename ' finish in ', num2str(toc), ' s']);