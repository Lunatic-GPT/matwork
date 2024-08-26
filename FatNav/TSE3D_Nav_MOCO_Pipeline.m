function TSE3D_Nav_MOCO_Pipeline(fname)

prefix=strtok2(fname,'.');

%SiemensRawData2mat(fname);
if 0% strcmp(rawDataVersion(fname),'vdve')
    a=mapVBVD(fname);
    maxLin=a{2}.hdr.MeasYaps.sWipMemBlock.alFree{64}-1;
    maxPar=a{2}.hdr.MeasYaps.sWipMemBlock.alFree{63}-1;
    
    mid=strtok_no(fname,'_',2);
    fmat=[mid,'_FatNav.mat'];

   if ~exist(fmat,'file')
    do_recon_FatNav_VE11(filename_append(fname(1:end-4),'_FatNav.mat'),maxLin,maxPar);
   end
    

  
    oprefix=strtok(fmat,'.');
   if ~exist([oprefix,'_vs.nii.gz'],'file')
    mat2nii_FatNav(fmat,fname);
    cmd = sprintf('3dvolreg -base %d -prefix %s_vs.nii.gz -dfile motion_%s.1D %s.nii',round(maxPar/2),oprefix,mid,oprefix);
    unix(cmd);  
  end

    cmd=sprintf('sbatch --mem=80G --ntasks=1 --job-name=a%s --time=64:00:00 mlb "run_grappa_and_moco_VE11(''%s'')"',mid,fname);
    unix(cmd);

    
else

 if ~exist([prefix,'_FatNav_recon.mat'],'file') 
  do_recon_FatNav([prefix,'_FatNav.mat'],[],[],false);
 end   
   if ~exist([prefix,'_vs.nii.gz'],'file')
    cmd = sprintf('3dvolreg -base 122 -prefix %s_vs.nii.gz -dfile motion_%s.1D %s_FatNav_recon.nii',prefix,prefix,prefix);
    unix(cmd);
   end

 
    prefix_out=[prefix,'_recon'];
    prefix_out2=[prefix,'_recon_moco'];
    
    fpro=name4pat(fullfile(prefix,'*.pro'),1);
    kdata=recon_tse_vfl_Siemens(fname,prefix_out,fpro);

    dfile=sprintf('motion_%s.1D',prefix);
   
    moco_tse_vfl_afterGrappa(kdata,fpro,dfile,8,prefix_out2);

end



