function mat2niigz(fname,var,prefix,zip,parent_nii)
% parent_nii can be nifti or a orient structure
if ~exist('zip','var')
    zip=true;
end

if ~isa(fname,'char')
    d=fname;
else
    
    if ~exist('var','var')
        d=ri_d1(fname);
    else
        d=ri_d1(fname,[],[],var);
    end
end

if isa(d,'logical') || isa(d,'uint16')
    d=int16(d);
end

if ~exist('prefix','var')||isempty(prefix)
  prefix=strtok2(fname,'.');
end

if ~exist('parent_nii','var') 
    
    nii=make_nii(d);
    try
        orient=get_orient(fname);
        nii=nii_from_orient(nii,orient);
    catch      
        nii.untouch=1;
        nii.hdr.hist.magic='n+1';
    end
elseif ~isa(parent_nii,'char')
    nii=make_nii(d); 
    nii=nii_from_orient(nii,parent_nii);
else    
    nii=load_untouch_niigz(parent_nii);   
    nii.img=d;
end

  if zip
    save_untouch_niigz(nii,prefix);  
  else
    save_untouch_nii(nii,prefix); 
  end
  
    
    
