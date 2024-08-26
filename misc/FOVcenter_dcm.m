function FOVcenter_dcm(dname)
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

for i=1:length(dir_str)
    h=dicominfo(fullfile(dname,dir_str(i).name));
    
    if i==1
     pos=h.ImagePositionPatient;
     pix=h.PixelSpacing;
     mtx=h.AcquisitionMatrix;
     mtx(mtx==0)=[];
     disp(pos(1:2));
     disp(double(mtx).*pix);
     
    end
    tr(i)=h.AcquisitionNumber;
    pos_tmp=h.SliceLocation;
    if ~any(pos_tmp==pos)
      pos(end+1)=pos_tmp;
    %else
     %   break;
    end
           
end

pos2=unique(pos);
ns=length(pos2);


fprintf('offset = %f; nslice = %d; thickness = %f\n',mean(pos2),ns,abs(pos2(2)-pos2(1)));
    disp([mfilename ' finish in ', num2str(toc), ' s']);
    
    