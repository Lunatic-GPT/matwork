function sdt2afni3d(prefix,fid_dir)
% sdt2afni(prefix[,fid_dir])
% if prefix is not used, then it will be the same as the fid_dir
% first dim: ro; second dim: pe; third dim: slice;
% 1/22/2011: Wrote it. Tested for single-slice 2D data. XZ
% nophase: do not save phase image.default true.


if ~exist('fid_dir','var')
  fid_dir=[prefix,'.fid'];
end

img=rdSdt(prefix);
%refpos = readPar(fid_dir,'ref_pos')+1;
%img(:,:,:,refpos)=[];

lro = readPar(fid_dir,'lro');%in cm
lpe = readPar(fid_dir,'lpe');%in cm
lpe2 = readPar(fid_dir,'lpe2');%in cm

tr = readPar(fid_dir,'tr');
%img = slice_reorder(img,pss);    

    pro=readPar(fid_dir,'pro'); %in cm
    ppe=readPar(fid_dir,'ppe'); %in cm
ppe2=readPar(fid_dir,'ppe2'); %in cm

    sz=size(img);
    delta = [lro*10/sz(1),lpe*10/sz(2),lpe2*10/sz(3)];  % thickness in mm
    
    orig_delta(1,:) = [pro,ppe,ppe2]*10-delta.*(sz(1:3)/2-[0.5,0.5,0.5]);
    orig_delta(2,:) = delta;
  %   [img,orig_delta] = reorient_data(img,orig_delta,orient(2:end-1));
  %   should not need to reorient if data is from sdt file.  6/4/2012
     
    orig_delta(:,[1,3])=-orig_delta(:,[1,3]);   
    info.ORIGIN = orig_delta(1,:);
    info.DELTA = orig_delta(2,:);
    info.TAXIS_FLOATS = [0,tr,0,0,0];
    info.ORIENT_SPECIFIC = [1 3 5]; %L2R,A2P,S2I
    write_afni(img,info,prefix);
   
if ndims(img)==4
    fprintf('std2afni finished. Image dimensions:%d*%d*%d*%d\n',size(img));
else
    fprintf('std2afni finished. Image dimensions:%d*%d*%d\n',size(img));
end


