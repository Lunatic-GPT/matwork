function gui_FatNav_callback(params,sbutton)

if ~exist('sbutton','var')
    sbutton=get(gco,'String');
end

sid_orig=get(params,'sid');

if ~isempty(sid_orig)
    sid=textscan(sid_orig,'%s','delimiter',',');
    params=set(params,'sid','');
    
    for i=1:length(sid{1})
        
        cd(sid{1}{i});
        gui_FatNav_callback(params,sbutton);
        disp(['inside ',pwd]);
        
        cd('..');
    end
    
    disp(pwd);
    set(params,'sid',sid_orig);
    
    return;
end

if ~isempty(strfind(sbutton,'Extract Data'))
    
    f_raw=get(params,'raw file');
    
    fsb_all=dir(f_raw);
   
    dname=filename(pwd);
          
    newdname=fullfile('..',sprintf('%s_mat',dname));
    if ~exist(newdname,'dir')
        mkdir(newdname);
    end
    
    for j=1:length(fsb_all)
            prefix=strtok(fsb_all(j).name,'.');
            
            %if ~exist(fullfile(newdname,[prefix,'_FatNav.mat']),'file') || ~exist(fullfile(newdname,[prefix,'.mat']),'file')
                
             if ~exist(fullfile(newdname,[prefix,'.mat']),'file')
                
                [Data,PMUTimeStamp,freePara,navData]= extractNavData(fsb_all(j).name,Inf);
                
                movefile([prefix,'.mat'],newdname);
                movefile(prefix,newdname);
                if ~isempty(navData)
                    movefile([prefix,'_FatNav.mat'],newdname);
                end
            end
    end
    
  %  gui_FatNav_callback(params,'gen bsub');
  %  gui_FatNav_callback(params,'push to killdevil/longleaf and run');
end

if ~isempty(strfind(sbutton,'gen bsub'))
    
    f_raw=get(params,'raw file');
    
    fsb_all=dir(f_raw);
    step=get(params,'images per process');
    fatnavrecon_name='FatNavRecon.sh';
    fid_fatnav=fopen(fatnavrecon_name,'w');
    dname=filename(pwd);
    
         
    newdname=fullfile('..',sprintf('%s_mat',dname));
    if ~exist(newdname,'dir')
        mkdir(newdname);
    end
    
    for j=1:length(fsb_all)
       
        
        prefix=strtok(fsb_all(j).name,'.');
      
        f_fatnav=fullfile(newdname,[prefix,'_FatNav.mat']);
        if exist(f_fatnav,'file')
            
            load(f_fatnav,'RepetitionNav');
        pos_sag = readsPar(fullfile(newdname,prefix,[prefix,'.pro']),'sCuboid.sPosition.dSag');
        pos_cor= readsPar(fullfile(newdname,prefix,[prefix,'.pro']),'sCuboid.sPosition.dCor');
        pos_tra=readsPar(fullfile(newdname,prefix,[prefix,'.pro']),'sCuboid.sPosition.dTra');
        if isempty(pos_sag)
            pos_sag=0;
        end
        if isempty(pos_cor)
            pos_cor=0;
        end
        if isempty(pos_tra)
            pos_tra=0;
        end
        
        
        
            center=[pos_tra,pos_cor,pos_sag];
            mid=strtok_no(prefix,'_',2);
            % load([prefix,'_FatNav.mat'],'RepetitionNav');
            
            lend=max(RepetitionNav)+1;
            ind=1:step:lend;
            ind2=ind+step-1;
            if ind2(end)>lend
                ind2(end)=lend;
            end
       
        if 0 % kill devil decomissioned
            
            for i=1:length(ind)
                fprintf(fid_fatnav,'bsub -M 16 -o logFatNav%s.%%J matbgk "do_recon_FatNav(''%s_FatNav.mat'',%d:%d)" FatNav%s_%d.log\n',mid,prefix,ind(i),ind2(i),mid,i);
                
            end
            fprintf(fid_fatnav,'bsub -M 16 -o logvolreg%s.%%J matbgk "do_volreg_afterFatNav(''%s'',1:%d:%d,[%f,%f,%f])" volreg%s.log\n',mid,mid,step,lend,center,mid);
            
        else
            for i=1:length(ind)
                fprintf(fid_fatnav,'sbatch --mem=16000 --time=12:00:00 --job-name=FatNav%s matbgk "do_recon_FatNav(''%s_FatNav.mat'',%d:%d)" FatNav%s_%d.log\n',mid,prefix,ind(i),ind2(i),mid,i);
                
            end
            fprintf(fid_fatnav,'sbatch --mem=16000 --time=12:00:00 --job-name=volreg%s matbgk "do_volreg_afterFatNav(''%s'',1:%d:%d,[%f,%f,%f])" volreg%s.log\n',mid,mid,step,lend,center,mid);
            
        end
%             
%                    fprintf(fid,'#!/bin/bash\n');
%         
%                fprintf(fid,'#SBATCH --job-name=_%s\n',mid);
%         fprintf(fid,'#SBATCH --ntasks=1\n');
%         fprintf(fid,'#SBATCH --time=12:00:00\n');
%              fprintf(fid,'#SBATCH --mem=16000\n');

        
            if ~isempty(strfind(fsb_all(j).name,'tse_vfl'))
                raw_file=get(params,'raw file');
                dfile=get(params,'dfile');
                
                params=set(params,'raw file',[prefix,'.dat']);
                params=set(params,'dfile',sprintf('Motion_%s.1D',mid));
                gui_FatNav_callback(params,'recon T2 w/mc bsub');
                params=set(params,'raw file',raw_file);
                params=set(params,'dfile',dfile);               
            end
        end
        
       
    end
    fclose(fid_fatnav);
   
    copyfile(fatnavrecon_name,newdname);
end

if ~isempty(strfind(sbutton,'push to longleaf'))
   
    
    dname=pwd;
    sid=filename(dname);
    cmd='winscp.exe /command "open longleaf" "option batch on" "cd /pine/scr/z/o/zongx"';
    
    if ~dir_exist(sid,true)
    cmd2= sprintf('%s "mkdir %s"',cmd,sid);
    fprintf('"mkdir %s"\n',sid);
    fail_free(cmd2);
    end
    
    cmd=sprintf('%s "cd %s"',cmd,sid);
    cmd=sprintf('%s "lcd %s/../%s_mat"',cmd,pwd,sid);
    
    fsb_all=dir(sprintf('../%s_mat/*_FatNav.mat',sid));
 
    for j=1:length(fsb_all)
        if ~dir_exist([sid,'/',fsb_all(j).name],true)
          cmd2=sprintf('%s "put %s"',cmd,fsb_all(j).name);   
          fprintf('"put %s"\n',fsb_all(j).name); 
          fail_free(cmd2);
        end
        
        cmd2=sprintf('%s "put %s"',cmd,fsb_all(j).name(1:end-11));
        fprintf('"put %s"\n',fsb_all(j).name(1:end-11));
        fail_free(cmd2);
    end
    
    cmd2=sprintf('%s "put FatNavRecon.sh"',cmd);
    cmd2=sprintf('%s "exit"',cmd2);
    fprintf('"put FatNavRecon.sh"\n');
    fail_free(cmd2);
    
    % upload the protocol file
    fsb_all=dir(sprintf('../%s_mat/*_tse_vfl*.mat',sid));
    for j=1:length(fsb_all)     
       if ~strcmp(fsb_all(j).name(end-10:end),'_FatNav.mat') 
           mid=strtok_no(fsb_all(j).name,'_',2);
           cmd2=sprintf('%s "put recon_script_%s.sh"',cmd,mid);
           fprintf('"put recon_script_%s.sh"\n',mid);
           cmd2=sprintf('%s "exit"',cmd2);
           fail_free(cmd2);
           cmd2=sprintf('%s "put %s" "exit"',cmd,fsb_all(j).name(1:end-4));
           fprintf('"put %s"\n',fsb_all(j).name(1:end-4));
           fail_free(cmd2);
       end     
    end
    % upload the data
     for j=1:length(fsb_all)     
       if ~strcmp(fsb_all(j).name(end-10:end),'_FatNav.mat') 
           
           if ~dir_exist([sid,'/',fsb_all(j).name],false)
               
               cmd2 = sprintf('%s "put %s" "exit"',cmd,fsb_all(j).name);
               fprintf('"put %s"\n',fsb_all(j).name);
               fail_free(cmd2);
           end
  
       end     
     end
  
end
 if ~isempty(strfind(sbutton,'run'))
     
        %% push to longleaf
     
    sid=filename(pwd);
    cmd='winscp.exe /command ';
    cmd=[cmd,'"open longleaf" "option batch on" "cd /pine/scr/z/o/zongx"'];
   
    
    cmd=sprintf('%s "cd %s"',cmd,sid);
    
 cmd=sprintf('%s "call source FatNavRecon.sh"',cmd);
 
 fail_free(cmd);
 
    fsb_all=dir(sprintf('../%s_mat/*_FatNav.mat',sid));
  for j=1:length(fsb_all)     
       if ~strcmp(fsb_all(j).name(end-10:end),'_FatNav.mat') 
                     
           mid=strtok_no(fsb_all(j).name,'_',2);
           fprintf('"call sbatch recon_script_%s.sh"',mid);
           cmd2=sprintf('%s "call sbatch recon_script_%s.sh"',cmd,mid);
           fail_free(cmd2);
       end     
  end
     
end

if ~isempty(strfind(sbutton,'get output'))
    
    
    sid=filename(pwd);
    dir_str=dir('*_tse_vfl*.dat');
    dest=fullfile(todb,'PVS_R21_Temp',sid);
    mkdir(dest);
 % get motion files   
    for i=1:length(dir_str)
        prefix=strtok(dir_str(i).name,'.');
        mid=strtok_no(prefix,'_',2);
        
        dname='/nas/longleaf/home/zongx/MotionFiles';
        while ~dir_exist([dname,'/',sid],false)
            
            fprintf('waiting for folder %s\n',sid);
            pause(3600);
        end
        
        fname = sprintf('Motion_%s.1D',mid);
            
        while ~dir_exist([dname,'/',sid,'/',fname],false)
           
            fprintf('waiting for file %s\n',fname);
            
            pause(3600);
        end
        winscp_get(0,[dname,'/',sid],fname);
        
        movefile(fname,dest);
        
    end
 % get recon files
 
 
 for i=1:length(dir_str)
     prefix=strtok(dir_str(i).name,'.');
     
     dname='/nas/longleaf/home/zongx/pine';

     while ~dir_exist(sid,false)        
         fprintf('waiting for folder %s\n',sid);
          pause(3600);
     end
     
     fname = sprintf('%s_recon_MotionCorr_iter5.mat',prefix);
     
     while ~dir_exist([sid,'/',fname],false)
         fprintf('waiting for file %s\n',fname);
         
         pause(3600);
     end
     winscp_get(0,[dname,'/',sid],fname);
     movefile(fname,dest);
 end
 
end

if strcmp(sbutton,'recon')
    
fname=get(params,'raw file');
fname=filename(fname);
prefix=strtok(fname,'.');

    s=get(params,'Repititions: start stop');
    do_recon_FatNav([prefix,'_FatNav.mat'],s(1):s(2));
elseif strcmp(sbutton,'mat2afni')
    

fname=get(params,'raw file');
fname=filename(fname);
prefix=strtok(fname,'.');    
    mid=strtok_no(prefix,'_',2);
    s=get(params,'Repititions: start');
    mat2afni_FatNav([mid,'_FatNav'],s); %the third parameter for center is missing. but the center value does not matter since 3dvolreg is wrt the FOV center
elseif strcmp(sbutton,'volreg')
    
fname=get(params,'raw file');
fname=filename(fname);
prefix=strtok(fname,'.');
    mid=strtok_no(prefix,'_',2);
    cmd=sprintf('3dvolreg -dfile Motion_%s.1D -base 0 -prefix %s_FatNav_volreg %s_FatNav+orig',mid,mid,mid);
    disp(cmd);
    unix(cmd);
elseif strcmp(sbutton,'plot motion')
    dfile=get_fpattern(params,'dfile');
    dfile=str2cell(dfile);
    for i=1:length(dfile)
        if ~strcmp(dfile{i}(end-2:end),'mat')
          plot_motion_afni(dfile{i});
        else
          plot_motion_afni_newBase(dfile{i});      
        end
        ttl=filename(dfile{i});
        ttl=strtok(ttl,'.');
        ttl=strrep(ttl,'_',' ');
        title(ttl);
    end
elseif strcmp(sbutton,'ESPIRiT bsub')
    lBaseResolution = readsPar(fullfile(prefix,[prefix,'.pro']),'lBaseResolution');
    lPhaseEncodingLines =readsPar(fullfile(prefix,[prefix,'.pro']),'lPhaseEncodingLines');
    lPartitions = readsPar(fullfile(prefix,[prefix,'.pro']),'lPartitions');
    lImagesPerSlab = readsPar(fullfile(prefix,[prefix,'.pro']),'lImagesPerSlab');
    lMatrix=[ lBaseResolution, lPhaseEncodingLines, lImagesPerSlab,lPartitions ];
    
    fid=fopen('ESPIRiT_script.sh','w');
    nmaps=get(params,'maps');
    ind=get(params,'iro indices (start stop step)');
    ind=linspace(ind(1)-1,ind(2),round((ind(2)-ind(1)+1)/ind(3))+1);
    ind=round(ind);
    for i=1:length(ind)-1
        fprintf(fid,'bsub -M 26 -o logESPIRiT%i.%%J matbgk "recon_ESPIRiT3D(''%s.mat'',%d:%d,%d,[%d,%d,%d,%d])" logESPIRiT_mat%d\n',i,prefix,ind(i)+1,ind(i+1),nmaps,lMatrix,i);
    end
    
    str=num2str(ind+1);
    str=strrep(str,'  ',',');
    str=strrep(str,',,',',');
    
    fprintf(fid,'bsub -M 16 -o logESPIRiTCombine.%%J matbgk "do_combine_afterESPIRiT3D(''%s'',%d,''%s'')" logESPIRiTCombine_mat\n',prefix,nmaps,str);
    fclose(fid);
    
    
elseif strcmp(sbutton,'recon T2 w/mc bsub')
    
    
fname=get(params,'raw file');
fname=filename(fname);
prefix=strtok(fname,'.');


    dfile=get(params,'dfile');

    dname=filename(pwd);
          
    newdname=fullfile('..',sprintf('%s_mat',dname));
    
    fov(1)=readsPar(fullfile(newdname,prefix,[prefix,'.pro']),'asSlice[0].dReadoutFOV');
    fov(3)=readsPar(fullfile(newdname,prefix,[prefix,'.pro']),'asSlice[0].dThickness');
    fov(2)=readsPar(fullfile(newdname,prefix,[prefix,'.pro']),'asSlice[0].dPhaseFOV');
    
    
    lBaseResolution = readsPar(fullfile(newdname,prefix,[prefix,'.pro']),'lBaseResolution');
    lPhaseEncodingLines =readsPar(fullfile(newdname,prefix,[prefix,'.pro']),'lPhaseEncodingLines');
    lPartitions = readsPar(fullfile(newdname,prefix,[prefix,'.pro']),'lPartitions');
    lImagesPerSlab = readsPar(fullfile(newdname,prefix,[prefix,'.pro']),'lImagesPerSlab');
    lAccelFactPE =readsPar(fullfile(newdname,prefix,[prefix,'.pro']),'lAccelFactPE');
    lMatrix=[ lBaseResolution, lPhaseEncodingLines, lImagesPerSlab,lPartitions ];
    nmapjobs = get(params,'number of map calc jobs');
    nnufftjobs = get(params,'number of nufft blocks');
    
    mid=strtok_no(prefix,'_',2);
    fid=fopen(sprintf('recon_script_%s.sh',mid),'w');
    killDevil = false;
    if killDevil  % not stable on killDevil
        fprintf(fid,'bsub -M 16 -o logmap.%%J matbgk "prep_bsub_ESPIRiT3D_maps(''%s.mat'',''%s'',[%4.1f,%4.1f,%4.1f],[%d,%d,%d,%d],%d)" logmap_mat\n',prefix,dfile,fov,lMatrix,nmapjobs);
        %    recon_ESPIRiT3D_FatNav(fname,dfile,fov,lMatrix,njobs,file_res0)
        % failed at 80 MG
        fprintf(fid,'bsub -M 100 -o logrecon.%%J matbgk "recon_ESPIRiT3D_FatNav(''%s.mat'',''%s'',[%4.1f,%4.1f,%4.1f],[%d,%d,%d,%d],%d,%d)" logrecon_mat\n',prefix,dfile,fov,lMatrix,nnufftjobs,nmapjobs);
    else
        fprintf(fid,'#!/bin/bash\n');
        
        fprintf(fid,'#SBATCH --job-name=recon_prep_%s\n',mid);
        fprintf(fid,'#SBATCH --ntasks=1\n');
        fprintf(fid,'#SBATCH --time=24:00:00\n');
        fprintf(fid,'#SBATCH --mem=64000\n');
     %   fprintf(fid,'#SBATCH --partition=bigmem\n');
     %   fprintf(fid,'#SBATCH --qos bigmem_access\n');
               
   %  dfile=['/nas02/home/z/o/zongx/MotionFiles/',dname,'/',dfile];
   dfile=['~/MotionFiles/',dname,'/',dfile];
        if lAccelFactPE>1
        fprintf(fid,'matbgk "recon_ESPIRiT3D_FatNav(''%s.mat'',''%s'',[%4.1f,%4.1f,%4.1f],[%d,%d,%d,%d],%d,%d,'''',false)" logrecon_prep_%s\n',prefix,dfile,fov,lMatrix,nnufftjobs,nmapjobs,mid);
        else % when there is no undersampling
        fprintf(fid,'matbgk "recon_nufft_FatNav(''%s.mat'',''%s'',[%4.1f,%4.1f,%4.1f],[%d,%d,%d,%d])" logrecon_%s\n',prefix,dfile,fov,lMatrix,mid);
        end
    end
    fclose(fid);
    
    movefile(sprintf('recon_script_%s.sh',mid),newdname);
end

disp('gui_FatNav_callback done!');

function results=fail_free(cmd)
    status=1;
    n=0;
    while status~=0
       fprintf('...');
       disp(cmd);
       if n>=10
       disp('close winscp before running this; otherwise error');
       end
       n=n+1;
       [status,results]=system(cmd);
      % fprintf('+\n');
    end
    
    
        
        