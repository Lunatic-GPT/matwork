function sdt2afni(prefix,fid_dir)
% sdt2afni(prefix[,fid_dir])
% if prefix is not used, then it will be the same as the fid_dir
% first dim: ro; second dim: pe; third dim: slice;
% 1/22/2011: Wrote it. Tested for single-slice 2D data. XZ
% nophase: do not save phase image.default true.


if ~exist('fid_dir','var')
  fid_dir=[prefix,'.fid'];
end

img=rdSdt(prefix);
ns =  readPar(fid_dir,'ns');
%refpos = readPar(fid_dir,'ref_pos')+1;
%img(:,:,:,refpos)=[];

lro = readPar(fid_dir,'lro');%in cm
lpe = readPar(fid_dir,'lpe');%in cm
thk = readPar(fid_dir,'thk');%in mm
tr = readPar(fid_dir,'tr');
pss=readPar(fid_dir,'pss'); %in cm
orient=readPar(fid_dir,'orient');
%img = slice_reorder(img,pss);    

    pro=readPar(fid_dir,'pro'); %in cm
    ppe=readPar(fid_dir,'ppe'); %in cm
    sz=size(img);
    delta = [lro*10/sz(1),lpe*10/sz(2),thk];  % thickness in mm
    if length(pss)>1
        pss_sort = sort(pss);
        delta(3) = (pss_sort(2)-pss_sort(1))*10;  %pss in cm 
    end
    
    orig_delta(1,:) = [pro,ppe,min(pss)]*10-delta.*([sz(1:2)/2,0]-[0.5,0.5,0]);  
    
    %varian convention:
    %positive pro: axial P; axial90 R; coronal I; coronal 90 R; sag I; sag90 P;
    %positive ppe: axial L; axial90 P; coronal L; coronal90 I; sag A; sag90 I;
    %
    % 
    
    orig_delta(2,:) = delta;                                                     %    
  %   [img,orig_delta] = reorient_data(img,orig_delta,orient(2:end-1));
  %   should not need to reorient if data is from sdt file.  6/4/2012
     
    orig_delta(:,[1,3])=-orig_delta(:,[1,3]);   
    info.ORIGIN = orig_delta(1,:);
    info.DELTA = orig_delta(2,:);
    info.TAXIS_FLOATS = [0,tr,0,0,0];
    info.ORIENT_SPECIFIC = [1 3 5]; %L2R,A2P,S2I
    write_afni(img,prefix,info);
   
if ndims(img)==4
    fprintf('std2afni finished. Image dimensions:%d*%d*%d*%d\n',size(img));
elseif ndims(img)==3
    fprintf('std2afni finished. Image dimensions:%d*%d*%d\n',size(img));
else
    fprintf('std2afni finished. Image dimensions:%d*%d\n',size(img));
end


