function [d,voxsize]=ri_mat_dcm(fname,ns)


if exist([fname,'.mat'],'file')
    d=ri([fname,'.mat'],'','','d');
    
    voxsize=ri([fname,'.mat'],'','','voxsize');
else
  %  d=dcm2mat_siemens(swiname);
  if exist('ns','var')
      d=ri(fname,ns);
  else
      d=ri(fname);
  end
    voxsize=dcmDimCenter(fname);
    save([fname,'.mat'], 'd','voxsize');
end
