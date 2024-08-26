function dcm2enhanced(dname)
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
    brikdata(:,:,1,i) = dicomread(fullfile(dname,dir_str(i).name));
        
   
end

pos2=unique(pos);
brikdata=reshape(brikdata,[size(brikdata,1),size(brikdata,2),length(pos2),size(brikdata,4)/length(pos2)]);

mkdir(sprintf([dname,'_enhanced']));
cd(sprintf([dname,'_enhanced']));


dicomwrite(uint16(brikdata),sprintf('%s.dcm',dname),dinfo(1), 'ObjectType', 'Secondary Capture Image Storage' );

cd('..');
%% save the average time series for each trial type.    

    disp([mfilename ' finish in ', num2str(toc), ' s']);