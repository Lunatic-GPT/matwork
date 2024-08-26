function [d,voxsize]=ri_mat_dcm_ana(fname,ns)

[d,prefix]=fileparts(fname);


prefix=fullfile(d,prefix);

if exist([prefix,'.mat'],'file')
    d=ri([prefix,'.mat'],'','','d');
    
    voxsize=ri([prefix,'.mat'],'','','voxsize');
elseif exist([prefix,'.hdr'],'file')
    
    b=load_nii([prefix,'.hdr']);
    
    voxsize=b.hdr.dime.pixdim(2:4);
    d=b.img;
    
  %  d=dcm2mat_siemens(swiname);
else
    
    d=ri(fname,ns);
    voxsize=dcmDimCenter(fname);
    save([fname,'.mat'], 'd','voxsize');
end
