function gui_veinMoment_callback(params,sbutton,TE)

if ~exist('sbutton','var')
  sbutton=get(gco,'String');
end

if isa(params,'char')
    load(params);
end

if ~isempty(strfind(sbutton,'calc moment (Frangi)'))
    
     mag=get_fpattern(params,'mag dcm dir');
     ph=get_fpattern(params,'ph dcm dir');
     mask=get_fpattern(params,'mask name');
     B0=get_fpattern(params,'B0 (T)');
     if ~exist('TE','var')
       TE = readdPar(ph,'EchoTime');
     end
     vox_size=dcmDimCenter(ph);
     
     %%
   if  debug
     ph='TE2';
     mask='mask_wm.nii';
    mag='../MAG_IMAGES_0033/TE2';
    TE = 15;
    B0=7;
     vox_size=[0.4297    0.4297    0.400];
   end
      Calc_MagMoment_Frangi3D(ph,mag,mask,TE,B0,vox_size);
   
    
    %% debug
%     
%  ph='../PHA_IMAGES_0025/TE2';
%  mag='../MAG_IMAGES_0024/TE2';
%  TE = 14.3;
% B0=7;
% vox_size=[0.4,0.4,0.4];

elseif ~isempty(strfind(sbutton,'calc_moment_Neurite'))
    
    
    pathmask=get(params,'path2mask output');
  %  pathmask=filename(dir_trace);
   % maskfile=[pathmask,'.mat'];
    %V = get_path_direction(pathmask);  % already done 
    V = ri(pathmask,'','','V');
    mag=get_fpattern(params,'mag dcm dir');
    ph=get_fpattern(params,'ph dcm dir');
    %mask=get_fpattern(params,'mask name');
    B0=get_fpattern(params,'B0 (T)');
    if ~exist('TE','var')
        TE = readdPar(ph,'EchoTime');
    end
    if exist([ph,'.mat'],'file')
       vox_size=ri([ph,'.mat'],'','','voxsize');
    else
    vox_size=dcmDimCenter(ph);
    end
    
    Calc_MagMoment_Neurite(ph,mag,V,TE,B0,vox_size);
elseif strcmp(sbutton,'gen bsub script')
  
    save SaveParams4Bsub_veinMoment.mat params
    fid=fopen('calcMoment.sh','w');
    fprintf(fid,'bsub -M 4 -o logveinMoment.%%J matbgk "gui_veinMoment_callback(''SaveParams4Bsub_veinMoment.mat'',''calc_moment_Neurite'')" veinMoment.log\n');
          
    fclose(fid);
    

elseif ~isempty(strfind(sbutton,'show histogram'))
    
    load('res_all_center_neurite.mat','rad','dchi');
    dchi2=get(params,'dchi (ppm)');
    
    y=rad(rad>0)*2*sqrt(dchi/dchi2);
    
    figure;
    hist(y,30);
    set(gca,'FontSize',12);
    xlabel('Diameter (\mum)');
    ylabel('# of Voxels');
    xlim([10,50]);
    set(gca,'XTick',10:5:50);
    savetiffc Diameter_distr 
else
    
end
disp('gui_veinMomemt_callback finished');

function V = get_path_direction(maskfile)

 prefix=strtok(maskfile,'.');

    ThinningPathFind(maskfile,false);
    
    %% step 3b
    load(['ThPth_',maskfile]);  % need to fix connected paths afterwards
    %% check path length
%     for i=1:length(i_ind_path)
%         len=pathLength(i_ind_path{i},ind{i},voxsize,size(c));
%         fprintf('Path Length %d: %f\n',i,len);
%     end
%     
    
    pathCurved_nofix(['ThPth_',maskfile]);  %fixedpath
    
    prefix=strtok(maskfile,'.');
    fixGap(['ThPth_',prefix,'_Curv']);   %probv
    
   
    V=Path_Orientation_AllPathVox(['ThPth_',prefix,'_Curv_Dv.mat'],false);
    
    
 
    
    







disp('gui_spv_callback done');