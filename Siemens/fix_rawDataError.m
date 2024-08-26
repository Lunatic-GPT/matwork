flist={'meas_MID82_fl_fq_retroZ_mb_FatNav_FID24955.mat',...
    'meas_MID80_fl_fq_retroZ_sbRef_FatNav_FID24953.mat',...
    'meas_MID65_fl_fq_retroZ_ICA_FID24942.mat',...
    'meas_MID47_fl_fq_retroZ_FatNav_FID24924.mat',...
    'meas_MID40_fl_fq_retroZ_FatNav_FID24917.mat',...
    'meas_MID38_fl_fq_retroZ_FatNav_FID24915.mat'};

for i=1:length(flist)
    
    
    load(flist{i});
    
    n0=size(Data,2);
    n=floor(n0/32)*32;
    
    Data=Data(:,1:n);
    Line=Line(1:n);
    Slice=Slice(1:n);
    Partition=Partition(1:n);
    
    freePara = freePara(1:end-n0+n,:);
    PMUTimeStamp = PMUTimeStamp(1,1:end-n0+n);
    Set = Set(1,1:end-n0+n);
    
    name2=filename_append(flist{i},'2');
    save(name2,'Data','Line','Set','freePara','Partition','PMUTimeStamp','Slice');

end

    