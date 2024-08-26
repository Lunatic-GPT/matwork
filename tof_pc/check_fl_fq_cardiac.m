function hr=check_fl_fq_cardiac(fsb,nphase,f_fake,f_MON,maxRate,f_real,f_interp)


%nphase=get(params,'cardiac phases');
prefix=strtok(fsb,'.');
prefix_dir=fullfile(prefix,prefix);
if ~exist([prefix,'.mat'],'file')
    [Data,Line,Partition,Slice,Set,timeStamp,freePara,navData,LineNav,PartitionNav,SliceNav,RepetitionNav,dummyData]=readMeasDat(fsb,inf,0,true);
    szDummyData=size(dummyData);
    save([prefix,'.mat'],'Data','Line','Partition','Slice','Set','freePara','sizeDummyData');   
    if ~isempty(navData)
        save([prefix,'_FatNav.mat'],'navData','LineNav','PartitionNav','SliceNav','RepetitionNav');
    end    
else  
    load([prefix,'.mat'],'Line','Partition','Slice','Set','freePara'); 
    if ~exist('Line','var')% old format
        
        load([prefix,'.mat'],'lin','ushSet','freePara');
        Line=lin;
        Set=ushSet;
    end  
    dsb=zeros(1,32,length(Line)/32);  % dummy
end

nvenc=max(Set)+1;


%%


try
    phaseStab=readsPar([prefix_dir,'.pro'],'ucPhaseStabilize');
    
    if strcmp(phaseStab{1},'0x1')
        phaseStab=true;
    else
        phaseStab=false;
    end
catch
    phaseStab=false;
end

if  phaseStab
    
    dsb2=reorder_fl_fq_data_phaseStabOn(dsb,lin,prefix_dir,nvenc,do_retro,freePara(:,4),nphase,true);
else
    physio.f_fake=f_fake;
    physio.f_MON=f_MON;
    physio.f_real=f_real;
    physio.f_interp=f_interp;
    
    [dsb2,hr]=reorder_fl_fq_data(dsb,Line,prefix_dir,nvenc,freePara(:,4),nphase,physio,true,maxRate);
    
end
