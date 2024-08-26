function dcm_average(dname)
%dcm2analyze(dname,pattern,nslices,nt)
tic;
tr=[];
pos=[];

dir_str=dir([dname,'/*.dcm']);
if length(dir_str)==0
dir_str=dir([dname,'/*.IMA']);
end    

for i=1:length(dir_str)
    dinfo(i)=dicominfo(fullfile(dname,dir_str(i).name));

    tr(i)=dinfo(i).AcquisitionNumber;
    pos(i)=dinfo(i).SliceLocation;
    disp(i);
            img = dicomread(fullfile(dname,dir_str(i).name));
        
        
            
        brikdata(:,:,i) = img;

        
end

            vsize(1:2)=dinfo(1).PixelSpacing;
            vsize(3)=dinfo(1).SliceThickness;
            
pos2=unique(pos);
ns=length(pos2);

[tr2,ind]=sort(tr,'ascend');
brikdata=brikdata(:,:,ind);
dinfo=dinfo(ind);
d=reshape(brikdata,[size(brikdata,1),size(brikdata,2),ns,size(brikdata,3)/ns]);
d2=mean(d,4);
mkdir(sprintf([dname,'_mean']));
cd(sprintf([dname,'_mean']));

for i=1:size(d2,3)
 dicomwrite(uint16(d2(:,:,i)),sprintf('%s_%d.dcm',dname,i),dinfo(i));
end

cd('..');
%% save the average time series for each trial type.    

    disp([mfilename ' finish in ', num2str(toc), ' s']);