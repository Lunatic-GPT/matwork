function pos=dcm2ana_multichan(dname,nchan)
%dcm2analyze(dname,pattern,nslices,nt)


tic;

if ~exist('sl_first','var')
    sl_first=true;
end

tr=[];
img=[];
pos=[];

dir_str=dir([dname,'/*.dcm']);
if isempty(dir_str)
 dir_str=dir([dname,'/*.IMA']);
end    

zerobase=0;

for i=1:length(dir_str)
   if mod(i,10)==0
       disp(i);
   end
    
   if mod(i,nchan)==1 && ~exist('ns','var')
    h=dicominfo(fullfile(dname,dir_str(i).name));
   
        tmp=h.SliceLocation;
        if ~any(tmp==pos)
         pos(end+1)=tmp;
        else
            if nargout>0
                return;  %just get slice positions.
            end
            ns=length(pos);
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
        
        brikdata(:,:,i) = img;

        
end
if ~exist('ns','var')
    ns=length(pos);
end

d=reshape(brikdata,[size(brikdata,1),size(brikdata,2),nchan,ns,size(brikdata,3)/ns/nchan]);

d=permute(d,[1,2,4,3,5]);
d=reshape(d,[size(d,1),size(d,2),ns,size(d,4)*size(d,5)]);

d2=d;
[tmp,ind]=sort(pos);
d2=d(:,:,ind,:);

vsize(3)=tmp(2)-tmp(1);

writeanalyze(d2,dname,vsize,'int16');
     
%% save the average time series for each trial type.    

    disp([mfilename ' finish in ', num2str(toc), ' s']);