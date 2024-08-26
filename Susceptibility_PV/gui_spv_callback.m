function gui_spv_callback(params,sbutton)

if ~exist('sbutton','var')
    sbutton=get(gco,'String');
end

if isa(params,'char')  %for bsub
    load(params);
end

if ~isempty(strfind(sbutton,'Sep to '))
    % only for Siemens data
    nTE=sbutton(7:8);
    nTE = str2double(nTE);
    swi=get(params,'swi dir number');
    
%     
%     dirs{1}=str2cell(get_fpattern(params,'swi dcm dir'));
%     dirs{2}=str2cell(get_fpattern(params,'mag dcm dir'));
%     dirs{3}=str2cell(get_fpattern(params,'ph dcm dir'));
%     
dirs{2}=sprintf('MAG_IMAGES_%04d',swi-3);
dirs{1}=sprintf('SWI_IMAGES_%04d',swi);
dirs{3}=sprintf('PHA_IMAGES_%04d',swi-2);
dirs{4}=sprintf('MIP_IMAGES(SW)_%04d',swi-1);

    for i=1:4
            cd(dirs{i});
            d=dir2('*.IMA');
            if isempty(d)
                cd('..');
                continue;
            end
            for j=1:nTE
                mkdir(['TE',num2str(j)]);
                nfile=length(d)/nTE;
                for k=1:nfile
                    movefile(d((j-1)*nfile+k).name,['TE',num2str(j)]);
                end
                
            end
            
            cd('..');
            
    end
end

if ~isempty(strfind(sbutton,'mkdir swi_00#_ana'))
    swi=get(params,'swi dir number');
    mkdir(sprintf('swi_%04d_ana',swi));
end

if ~isempty(strfind(sbutton,'invert image'))
    
    swiname=get(params,'swi analyze file');
    
    d=load_untouch_niigz(swiname);
    
    prefix=strtok(swiname,'.');
    d.img=max(d.img(:))-d.img;
    save_untouch_niigz(d,[prefix,'_Inv']);
    disp('Invert Image Done!');
    
end
if ~isempty(strfind(sbutton,'path2mask enlarge'))
    
    
    swi_dir=filename(pwd);
    nswi=str2num(swi_dir(end-7:end-4));
    swiname=sprintf('../SWI_IMAGES_%04d/TE2.mat',nswi);
    
    dir_trace=get(params,'neurite path dir');
    
    [pathmask,name,suf]=fileparts(dir_trace);
    
    if isempty(pathmask)
        pathmask=dir_trace;
    else
        pathmask=[name,suf];
    end
    
   [tmp,maskfile]= add_neurite_traces(swiname,dir_trace,{'xy'},-2, pathmask,'');
    

    ThinningPathFind(maskfile,false);
    
    %% step 3b
    load(['ThPth_',maskfile]);  % need to fix connected paths afterwards
    %% check path length
    for i=1:length(i_ind_path)
        len=pathLength(i_ind_path{i},ind{i},voxsize,size(c));
        fprintf('Path Length %d: %f\n',i,len);
    end
    
    
    pathCurved_nofix(['ThPth_',maskfile]);  %fixedpath
    
    prefix=strtok(maskfile,'.');
    fixGap(['ThPth_',prefix,'_Curv']);   %probv
    
    %%
    
    tmp=load(['ThPth_',prefix,'_Curv_DV.mat']);
    V=calc_V(tmp.i_ind_path,tmp.ind,size(tmp.c));
    
    
    
    %%
    disp('determine mask path done!');
    
    
    disp('path2mask done!');
end

if ~isempty(strfind(sbutton,'path2mask')) && isempty(strfind(sbutton,'path2mask enlarge'))
    
    dir_trace=get(params,'neurite path dir');
    
    % pathmask=filename(dir_trace);
    
    
    [pathmask,name,suf]=fileparts(dir_trace);
    
    if isempty(pathmask)
        pathmask=dir_trace;
    else
        pathmask=[name,suf];
    end
    
    
    swi_dir=filename(pwd);
    nswi=str2num(swi_dir(end-7:end-4));
    swiname=sprintf('../PHA_IMAGES_%04d/TE2.mat',nswi-2);
   % add_neurite_traces(swiname,dir_trace,'*',{'xy'},[],pathmask,'');
    add_neurite_traces(swiname,dir_trace,{'xy'},[],pathmask,'');
    
    disp('path2mask done!');
end

if ~isempty(strfind(sbutton,'determine_mask_path'))
    
    maskfile=get(params,'mask file');
    
    ThinningPathFind(maskfile);
    
    %% step 3b
    load(['ThPth_',maskfile]);  % need to fix connected paths afterwards
    %% check path length
    for i=1:length(i_ind_path)
        len=pathLength(i_ind_path{i},ind{i},voxsize,size(c));
        fprintf('Path Length %d: %f\n',i,len);
    end
    
    
    pathCurved_nofix(['ThPth_',maskfile]);  %fixedpath
    
    prefix=strtok(maskfile,'.');
    fixGap(['ThPth_',prefix,'_Curv']);   %probv
    
    
    disp('determine mask path done!');
end
if ~isempty(strfind(sbutton,'rotate_data'))
    
    B0_dir=[0,0,1];
    crop_size=[64,64,32];  % crop size before interpolation.
    % output matrix size after rotate data.
    %lseg=3; % units mm.  path segment length.
    
    vox_size_new =get(params,'new vox size (mm)');
    
    lout = round(8/vox_size_new);
    
    lseg=get(params,'segment length (mm)');
    % lseg=3;
    
    
    maskfile=get(params,'mask file');
    
    
    [root,fname]=fileparts(maskfile);
    
    %maskfile=get(params,'mask path file');
    if ~isempty(root)
        cur_dir=cd(root);      
        prefix=strtok(fname,'.');
        
    else
        cur_dir=pwd;
        
        prefix=strtok(maskfile,'.');
    end
    
    dir4allvox=true;
    if lseg>0
      [p1,p2,an,resid_pathDir]=Path_Orientation(['ThPth_',prefix,'_Curv_Dv.mat'],lseg);
    else
     % [tmp,p1,p2,an,resid_pathDir]=Path_Orientation_AllPathVox(['ThPth_',prefix,'_Curv_Dv.mat']);
     
     [tmp,p1,p2,an,resid_pathDir]=Path_Orientation_AllPathVox(['ThPth_',prefix,'_Curv_Dv.mat'],dir4allvox);
    end
    
    
    csave = cell(1,length(an));
    for iTE=1:10
        try
            phdir=get_fpattern(params,sprintf('ph dcm dir (for TE%d)',iTE));
        catch
            break;
        end
        mgdir=get(params,sprintf('mag dcm dir (for TE%d)',iTE));
        
        if isempty(phdir)
            continue;
        end
        TE=readdPar(phdir,'EchoTime');
        phdcm=ri_mat_dcm(phdir,[]);
        magdcm=ri_mat_dcm(mgdir,[]);
        phrange=max(phdcm(:))-min(phdcm(:));
        vox_size=ri(['ThPth_',prefix,'_Curv_Dv.mat'],[],[],'voxsize');
        
        vox_size_new =get(params,'new vox size (mm)');
        interp=ceil(vox_size(1:3)*2/vox_size_new);
        data=double(magdcm).*exp(1i*double(phdcm)/double(phrange)*2*pi-1i*pi);
        
        lmax=maxlen_cell(an);
        
        mag=zeros(lout,lout,lmax,length(an));
        
        ph=zeros(lout,lout,lmax,length(an));
        deg=cell(1,length(an));
        
        for i=1:length(an)
            for j=1:length(an{i})
                %  prefix=sprintf('rotateDataOutput/path_%d_seg_%d',i,j);
                
                if (j<3 || j>length(an{i})-2)  && lseg==0%p1 p2 not well defined for end points
                    continue;
                end
                [mag_tmp,ph_tmp,deg{i}(j)]=rotate_data_Line2dim3_B2dim1(p1{i}(j,:),p2{i}(j,:),B0_dir,vox_size(1:3),data,crop_size,interp,vox_size_new);  %deg is the same as an
                c0=ceil(size(mag_tmp)/2+0.5);
                
                %c=centerFinder(mag_tmp(:,:,c0(3)).*exp(1i*ph_tmp(:,:,c0(3))),2,4);
                %      c1=centerFinder(mag_tmp(:,:,c0(3)),1,4,0);
                %      c2=centerFinder(mag_tmp(:,:,c0(3)),1,4,1);
                %      c3=centerFinder(mag_tmp(:,:,c0(3)).*exp(1i*ph_tmp(:,:,c0(3))),1,4,0);
                
                if iTE==1
                    csave{i}(j,:)=centerFinder(mag_tmp(:,:,c0(3)).*exp(1i*ph_tmp(:,:,c0(3))),1,4,1);
                end
                %  ctmp=centerFinderPhase(ph_tmp(:,:,c0(3)),4,i_rad_roi,i_rad_bg);
                
                c=csave{i}(j,:);
                mag(:,:,j,i) = mag_tmp(c(1)-lout/2:c(1)+lout/2-1,c(2)-lout/2:c(2)+lout/2-1,c0(3));
                ph(:,:,j,i) = ph_tmp(c(1)-lout/2:c(1)+lout/2-1,c(2)-lout/2:c(2)+lout/2-1,c0(3));
                
                %  mag(:,:,j,i) = mag_tmp(c0(1)-lout/2:c0(1)+lout/2-1,c0(2)-lout/2:c0(2)+lout/2-1,c0(3));
                %  ph(:,:,j,i) = ph_tmp(c0(1)-lout/2:c0(1)+lout/2-1,c0(2)-lout/2:c0(2)+lout/2-1,c0(3));
                
            end
        end
        
        save(['Rot_Ph_Mag_TE',num2str(iTE),'_',prefix,'_dir4allvox.mat'],'mag','ph','deg','resid_pathDir','TE','vox_size','vox_size_new');
        
    end
    cd(cur_dir);
    disp('rotate data done!');
end

if ~isempty(strfind(sbutton,'calc')) ||  ~isempty(strfind(sbutton,'calc_from_pattern'))
    
    
    te4moment = get(params,'TE for calc moment');
    
    te4s=get(params,'TE for calc susc');
    
    te4s=str2num(te4s);
    
    lseg=get(params,'segment length (mm)');
    rad_v = get(params,'vessel ROI radius (mm)');
    rad_mom = get(params,'Mom calc ring radii (mm)');
    rad_bg = get(params,'bg ring radii (mm)');
    
    rad_v = str2num(rad_v);
    
    
    B0=get(params,'B0 (T)');
    
    
    for iTE=1:max([te4s,te4moment])
        
        
   %    dname=get(params,sprintf('rotated data %d',iTE));
        mask_file=get(params,'mask file');
        fname=filename(mask_file);
        prefix=strtok(fname,'.');
        
        dname = ['Rot_Ph_Mag_TE',num2str(iTE),'_',prefix,'_dir4allvox.mat'];
    %  dname = ['Rot_Ph_Mag_TE',num2str(iTE),'_',prefix,'.mat'];
    
                debug=false;%true;
            
   
        vox_size=ri(dname,[],[],'vox_size');
        vox_size_new =ri(dname,[],[],'vox_size_new');
        
        TE(iTE)=ri(dname,'','','TE');
        mag(:,:,:,:,iTE)=ri(dname,'','','mag');
        ph(:,:,:,:,iTE)=ri(dname,'','','ph');
        if iTE==1
            deg=ri(dname,'','','deg');
            resid_pathDir=ri(dname,'','','resid_pathDir');
        end
        
    end
    
    cmplxd = mag.*exp(1i*ph);
    
    res=cell(1,length(deg));
    flag=cell(1,length(deg));
    cd_fit=cell(1,length(deg));
    cd_data=cell(1,length(deg));
    roi=cell(1,length(deg));
    
    resid_magMoment=cell(1,length(deg));
    resid=cell(1,length(deg));
    
    m_mom=mag(:,:,:,1,1)*0;
    m_v=mag(:,:,:,1,1)*0;
    m_bg=mag(:,:,:,1,1)*0;
    npath=length(deg);
    global plot_results
    plot_results=false;
    %  c=ceil((size(mag(:,:,1,1))+1)/2);
    
    %     c(:,5,1)=[39,43]';
    %     c(:,5,1)=[41,41]';
    %
    %      c(:,1:size(mag,3),1:npath)=repmat([41,41]',[1,size(mag,3),npath]);
    %
    c=repmat(ceil((size(mag(:,:,1,1))+1)/2)',[1,size(mag,3),npath]);
    %     c(:,5,1) = [19,22];
    %     c(:,5,1) = [19,22];
    
    i_rad_v=round(rad_v/vox_size_new);
    i_rad_mom=round(rad_mom/vox_size_new);
    i_rad_bg=round(rad_bg/vox_size_new);
    
    count = 0;
    vessel_seg=get(params,'vessel seg');
  
    if isempty(vessel_seg)
        vessels=1:npath;
        segments=1:max(cellLen(deg));
    else
        vessels=vessel_seg(1);
        segments=vessel_seg(2);
    end
    
    for i=vessels
        for j=segments
            
                if length(segments)>1 && (j<3 || j>length(deg{i})-2) && lseg==0 %p1 p2 not well defined for end points
                    continue;
                end

            m_vtmp=mask_circle(size(mag(:,:,1,1)),i_rad_v,c(:,j,i),1);
            m_momtmp = mask_circle(size(mag(:,:,1,1)),i_rad_mom(1),c(:,j,i),1);
            m_momtmp2 = mask_circle(size(mag(:,:,1,1)),i_rad_mom(2),c(:,j,i),1);
            
            m_bgtmp=mask_circle(size(mag(:,:,1,1)),i_rad_bg(1),c(:,j,i),1);
            m_bgtmp2=mask_circle(size(mag(:,:,1,1)),i_rad_bg(2),c(:,j,i),1);
            
            m_v(:,:,j,i)=m_vtmp;
            m_mom(:,:,j,i)=m_momtmp==0&m_momtmp2>0;
            m_bg(:,:,j,i)=m_bgtmp==0&m_bgtmp2>0;
            ph_bg=angle(mean_roi(cmplxd(:,:,j,i),m_bg(:,:,j,i)));
            
            if strcmp(sbutton,'calc')
                [res{i}(j,:),flag{i}(j),resid_magMoment{i}(j)]= calc_vessel_susceptibility(squeeze(cmplxd(:,:,j,i,:)).*exp(-1i*ph_bg),m_v(:,:,j,i),m_mom(:,:,j,i),m_bg(:,:,j,i),deg{i}(j),TE,B0,vox_size_new,te4moment,te4s);
                disp([res{i}(j,3),resid_magMoment{i}(j)]);
                
                
                mname='mask_vessels_calc.mat';
                save(mname,'m_v','m_mom','m_bg');
                
            else
                
                
                cmplxd2{i}(:,:,:,j) = squeeze(cmplxd(:,:,j,i,:)).*exp(-1i*ph_bg);
                
      
         %%       
                tic;
                disp([i,j]);
           
                if ~debug
                    lnew=[20,20];
                [res{i}(j,:),resid{i}(j),cd_fit{i}(:,:,:,j),cd_data{i}(:,:,:,j),roi_tmp] = calc_vessel_susceptibility_FromPattern( cmplxd2{i}(:,:,:,j),rad_v,deg{i}(j),TE,B0,vox_size(1:2),vox_size(1:2)./vox_size_new,te4s,[20,20]);
                else
                %load the results from old fit; old fit did not have roi
                %and cmplxd2; calculate here and save to a different file.
                tmp=load('ResultsPattern_septalVein_debug_7.6ms_15.0ms.mat');
                
                center=ceil(size(cmplxd(:,:,1,1,1)+1)/2)+tmp.res{i}(j,end-1:end)./vox_size_new;
                roiRad_i=round(mean(rad_v./vox_size_new));
                 
                roi_tmp=mask_circle(size(cmplxd(:,:,j,i,1)),roiRad_i,center,1);   
                
                res{i}(j,:)=tmp.res{i}(j,:);
                resid{i}(j)=tmp.resid{i}(j);
                cd_fit{i}(:,:,:,j)=tmp.cd_fit{i}(:,:,:,j);
                cd_data=tmp.cd_data;    
                end
                
                count = count+1;
                tmp=cell2array(deg);
                fprintf('Remaining time = %f\n',toc*(length(tmp)-count));
                
                cd_fit{i}(:,:,:,j)=setv_roi(cd_fit{i}(:,:,:,j),roi_tmp==0,0);
                cd_data{i}(:,:,:,j)=setv_roi(cd_data{i}(:,:,:,j),roi_tmp==0,0);
                roi{i}(:,:,j)=roi_tmp;
               
                if length(segments)==1
                    for icd=1:size(cd_fit{i},3) 
                       figure(103+icd);
                      compare_cd(cd_data{i}(:,:,icd,j), cd_fit{i}(:,:,icd,j));
                    
                      title(sprintf('D = %3.2f mm; chi = %3.2fppm',2*res{i}(j,1),res{i}(j,2)));
                      fprintf('Angle = %4.2f\n',deg{i}(j));
                    end
                    
                    
                end
            end
        end
    end
    
    

    resname=get_res_fname(params);
    if length(segments)>1
       save(resname,'res','vox_size_new','TE','deg','B0','flag','resid_magMoment','resid_pathDir','resid','cd_fit','cd_data','roi','cmplxd2');
    end
    
    disp('calc done');
end

if strcmp(sbutton,'show pattern')
     
    % PVS08; septal vein 1; segment 11 has nice matched pattern.
    
    
    resname=get_res_fname(params);
    load(resname);
    
    for i=1:2
        
       for j=1:size(cd_data{i},4)
          
           
          if any(vec(cd_fit{i}(:,:,:,j))~=0)
           disp([i,j]);
              figure(102);
              cd_data{i}=cd_data{i}*exp(-1i*res{i}(j,7));           
              cd_fit{i}=cd_fit{i}*exp(-1i*res{i}(j,7));
              imdata=cat(4,real(cd_data{i}(:,:,1,j)),imag(cd_data{i}(:,:,1,j)),real(cd_data{i}(:,:,2,j)),imag(cd_data{i}(:,:,2,j)));
              imfit=cat(4,real(cd_fit{i}(:,:,1,j)),imag(cd_fit{i}(:,:,1,j)),real(cd_fit{i}(:,:,2,j)),imag(cd_fit{i}(:,:,2,j)));
              imall=cat(4,imdata,imfit);
              imall=repmat2(imall,4);
              imshow4(imall,[],[2,4],1);
              
              figure(103);
              cd_data{i}=cd_data{i}*exp(-1i*res{i}(j,7));           
              cd_fit{i}=cd_fit{i}*exp(-1i*res{i}(j,7));
              
              mag=cat(4,abs(cd_data{i}(:,:,1,j)),abs(cd_data{i}(:,:,2,j)),abs(cd_fit{i}(:,:,1,j)),abs(cd_fit{i}(:,:,2,j)));
              ph=cat(4,angle(cd_data{i}(:,:,1,j)),angle(cd_data{i}(:,:,2,j)),angle(cd_fit{i}(:,:,1,j)),angle(cd_fit{i}(:,:,2,j)));
              
              ph=scale2n(ph,100,[-pi,pi]);
              mag=scale2n(mag,100,[0,max(mag(:))]);
              
              imall=cat(4,ph,mag);
              imall=repmat2(imall,4);
              imshow4(imall,[],[2,4],1);
              
              pause;
          end
      
          
       end
        
    end
    
end
if strcmp(sbutton,'gen bsub')
    
    save SaveParams4Bsub.mat params
    fid=fopen('calc_from_pattern.sh','w');
    fprintf(fid,'bsub -M 4 -o logspv.%%J matbgk "gui_spv_callback(''SaveParams4Bsub.mat'',''calc_from_pattern_rotate_data_determine_mask_path'')" spv.log\n');
    
    fclose(fid);
    
end
if strcmp(sbutton,'show results')
    
    resname=get(params,'saved results');
    
    load(resname);
    
    resid=cell2array(resid_magMoment,2);
    res2=cell2arrayMax(res,1);
    
    mom=res2(:,3)*1000;
    figure;hist(mom,0:1:15);
    xlim([0,15]);
    set(gca,'FontSize',12);
    ylabel('# of veins');
    xlabel('chi\timesr^2 (ppm mm^2)');
    
end

if strcmp(sbutton,'save params')
    
    p=params;
    root=fullfile(todb,'matwork','blat','EASY_GUI','savedParameters');
    save(fullfile(root,'gui_spv_params.mat'),'p');
    
end


disp('gui_spv_callback done');


function resname=get_res_fname(params)

        mask_file=get(params,'mask file');
        fname=filename(mask_file);
        prefix=strtok(fname,'.');
        
        
    te4moment = get(params,'TE for calc moment');
    
    te4s=get(params,'TE for calc susc');
    
    te4s=str2num(te4s);
 for iTE=1:max([te4s,te4moment])
        
        dname = ['Rot_Ph_Mag_TE',num2str(iTE),'_',prefix,'.mat'];   
        TE(iTE)=ri(dname,'','','TE');
      
 end
    
    
    
    TEstr=sprintf('_%5.1fms',TE(te4s));
    TEstr(TEstr==' ')=[];
    
  
        prefix=sprintf('ResultsPattern_%s_debug',strtok(mask_file,'.'));
    
        resname=sprintf('%s%s.mat',prefix,TEstr);
   
    function V=calc_V(i_ind_path,ind,sz)
    
    V=zeros([sz,3]);
    
    for j=1:length(i_ind_path)
        pos=ind2subb(sz,ind{j}(i_ind_path{j}));
        
        for i=1:length(i_ind_path{j})
        if i<3 || i>length(i_ind_path{j})-2
            continue;
        end
        
        pca_res=pca(pos(i-2:i+2,:));   
        V(pos(i,1),pos(i,2),pos(i,3),:) = pca_res(:,1);
        
        end
    end
    
    