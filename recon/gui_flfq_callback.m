function gui_flfq_callback(params,sbutton)

if ~exist('sbutton','var')
    sbutton=get(gco,'String');
end

%nf=get(params,'Files');


if strcmp(sbutton,'check pulse/gen bsub')
    
    %  for j=1:nf
    f_raw=get(params,sprintf('raw data'));
    %     if isempty(f_raw)
    %         continue;
    %     end
    nphase =get(params,'cardiac phases');
    f_fake=get(params,'fake peaks');
    f_MON=get(params,'MON peaks');   % the peak should not be one of the fake peaks. It will be removed from the fake peak if it is.
    f_real=get(params,'real peaks');
    fsb_all=dir(f_raw);
    maxRate=get(params,'max heart rate (per min)');
    
    f_interp=get(params,'interp n peaks');
    for i=1:length(fsb_all)
        fsb=fsb_all(i).name;
        hr=check_fl_fq_cardiac(fsb,nphase,f_fake,f_MON,maxRate,f_real,f_interp);
        mid=strtok_no(fsb,'_',2);
        save(['HR_',mid,'.mat'],'hr');
          
    end
  
    gui_flfq_callback(params,'gen recon bsub');
    % end
    
    
elseif strcmp(sbutton,'load physio')
    [files, dataDir] = uigetfile('physio*.mat','Load physio',pwd,'MultiSelect','off');
    
   tmp= load(fullfile(dataDir,files));
    
      params=set(params,'fake peaks',tmp.physio.f_fake);
      params=set(params,'MON peaks',tmp.physio.f_MON);   % the peak should not be one of the fake peaks. It will be removed from the fake peak if it is.
      params=set(params,'real peaks',tmp.physio.f_real);
      params=set(params,'max heart rate (per min)',tmp.physio.maxRate);
    
      params=set(params,'interp n peaks',tmp.physio.f_interp);
elseif strcmp(sbutton,'gen recon bsub')
    f_raw=get(params,'raw data');
    nphase =get(params,'cardiac phases');
  
   
    interp_factor=get(params,'interp');
    
    f_raw=filename(f_raw);
    mid=strtok_no(f_raw,'_',2);
  
    fphysio=['physio_',mid,'.mat'];
  
    fid=fopen(sprintf('bsub_%s.m',mid),'w');
    % if isempty(f_sbRef)
    % recon_fl_fq(fsb,interp_factor,nphase,physio,DownsampleFactor,do_scale)
    nphase =get(params,'cardiac phases');
    f_fake=get(params,'fake peaks');
    f_MON=get(params,'MON peaks');   % the peak should not be one of the fake peaks. It will be removed from the fake peak if it is.
    f_real=get(params,'real peaks');
    maxRate=get(params,'max heart rate (per min)');
    
    f_interp=get(params,'interp n peaks');
    physio.f_fake=f_fake;
    physio.f_real=f_real;
    physio.f_MON=f_MON;
    physio.f_interp = f_interp;
    physio.maxRate=maxRate;
    if exist(fphysio,'file')
        fphysio=inputdlg('physio name','physio name',1,{[fphysio,'a']});
        fphysio=fphysio{1};
    end
    save(fphysio,'physio');
    
    fprintf(fid,'recon_fl_fq(''%s'',[%s],%d,''%s'')',f_raw,num2str(interp_factor),nphase,fphysio);
    % end
    fclose(fid);
    
    save(['gui_params_gen_recon_bsub_',mid,'.mat']);
    
elseif strcmp(sbutton,'raw to mat')
    f_raw=get(params,sprintf('raw data'));
    
    fsb_all=dir(f_raw);
    
    
    dname=filename(pwd);
    newdname=fullfile('..',sprintf('%s_mat',dname));
    if ~exist(newdname,'dir')
        mkdir(newdname);
    end
    
    for i=1:length(fsb_all)
        
        fsb=fsb_all(i).name;
        prefix=strtok(fsb,'.');
        [dsb,lin,par,sl,ushSet,timeStamp,freePara,navData,LineNav,PartitionNav,SliceNav,RepetitionNav,dummyData]=readMeasDat(fsb,inf,0,true);
        sizeDummyData=size(dummyData);
        save([prefix,'.mat'],'dsb','lin','par','sl','ushSet','freePara','sizeDummyData');
        if ~isempty(navData)
            save([prefix,'_FatNav.mat'],'navData','LineNav','PartitionNav','SliceNav','RepetitionNav');
        end
        movefile([prefix,'.mat'],newdname);
        movefile(prefix,newdname);
        if ~isempty(navData)
            movefile([prefix,'_FatNav.mat'],newdname);
        end
        
    end
    
elseif strcmp(sbutton,'FatNav recon script')
    fname=get(params,sprintf('raw data'));
    fname=filename(fname);
    prefix=strtok(fname,'.');
    bsub_FatNav([prefix,'_FatNav.mat']);
    
elseif strcmp(sbutton,'recon')
    f_raw=get(params,'raw data');
    nphase =get(params,'cardiac phases');
    physio.f_fake=get(params,'fake peaks');
    physio.f_real=get(params,'real peaks');
    
    physio.f_MON=get(params,'MON peaks');
    interp_factor=get(params,'interp');
    f_all=dir(f_raw);
    f_sbRef=get(params,'SB ref');
    %    nonzero_fraction = get(params,'nonzero fraction');  % the fraction of k-space data to keep at the center for recon
    maxRate=get(params,'max heart rate (per min)');
    for i=1:length(f_all)
        %recon_fl_fq(fsb,interp_factor,nphase,physio,DownsampleFactor)
        
        if isempty(f_sbRef)
            recon_fl_fq(f_all(i).name,interp_factor,nphase,physio);
        else
            recon_fl_fq_mb_grappa(f_sbRef,f_all(i).name,interp_factor,nphase,f_fake,f_MON,maxRate,f_real); % fix later
        end
        
        
    end
    
elseif strcmp(sbutton,'gen mb recon bsub')
    maxRate=get(params,'max heart rate (per min)');
    f_raw=get(params,'raw data');
    nphase =get(params,'cardiac phases');
    f_fake=get(params,'fake peaks');
    f_MON=get(params,'MON peaks');
    interp_factor=get(params,'interp');
    
    f_sbRef=get(params,'SB ref');
    f_raw=filename(f_raw);
    mid=strtok_no(f_raw,'_',2);
    
    fid=fopen(sprintf('bsub_%s.sh',mid),'w');
    if ~isempty(f_sbRef)
        fprintf(fid,'bsub -M 8 -o logflfqmbRecon_%s.%%J matbgk "recon_fl_fq_mb_grappa(''%s'',''%s'',[%f,%f],%d,''%s'',''%s'',%d)" logflfqRecon_%s',mid,f_sbRef,f_raw,interp_factor,nphase,f_fake,f_MON,maxRate,mid);
    end
    fclose(fid);
    
    %elseif strcmp(sbutton,'Clear')
    
    
end

disp('gui_flfq_callback Done!');


