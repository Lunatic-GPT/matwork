function nifti_resample(fname,fout,voxsize,matsize,use_nn)

    d=load_untouch_niigz(fname);
    
    pixdim_old=d.hdr.dime.pixdim(2:4);
    dim_old=d.hdr.dime.dim(2:4);
  
    img2=reslice(double(d.img),eye(3),d.hdr.dime.pixdim(2:4),[0,0,0],matsize,voxsize,false,use_nn);
    d.hdr.dime.dim(2:4)=matsize;
    d.hdr.dime.pixdim(2:4)=voxsize;
    d.img=int16(img2);
    
    srow=[d.hdr.hist.srow_x;d.hdr.hist.srow_y;d.hdr.hist.srow_z];
    
    srow2=srow(:,1:3).*(voxsize./pixdim_old)';
    
    center=pos2center(srow,dim_old);
    
    srow2(:,4)=center2pos(srow2(:,1:3),center,matsize);
    d.hdr.hist.srow_x=srow2(1,:);
    d.hdr.hist.srow_y=srow2(2,:);
    d.hdr.hist.srow_z=srow2(3,:);
    d.hdr.hist.qoffset_x=srow2(1,4);
    d.hdr.hist.qoffset_y=srow2(2,4);
    d.hdr.hist.qoffset_z=srow2(3,4);
    
    
    
    save_untouch_niigz(d,fout);
    
    
    
function center=pos2center(srow,sz)

 
    ijk_c=floor(double(sz)/2);      
      
    center=srow(:,4)+srow(:,1:3)*ijk_c(:);
    
    
    
function pos=center2pos(srow_1_3,center,sz)

 
    ijk_c=floor(double(sz)/2);      
      pos=  center-srow_1_3*ijk_c(:);
    
