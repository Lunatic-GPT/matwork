function Data = readMeasDat_DTI(filename,slice)


prefix=strtok(filename,'.');
fid = fopen(filename,'rb');
nProtHeaderLen = double(fread(fid,1,'*int32'));

if exist([filename,'_protocol'],'file')
   fseek(fid, nProtHeaderLen - 4, 0);
else
  header=fread(fid,nProtHeaderLen-4,'*char');
  fid2=fopen([filename,'_protocol'],'w');

  for i=5:length(header)-12
   fprintf(fid2,'%c',header(i));
  end
  fclose(fid2);

  extract_protocol([filename,'_protocol']);
end

mtx(1)=readsPar([prefix,'.pro'],'lBaseResolution');
mtx(2)=readsPar([prefix,'.pro'],'lPhaseEncodingLines');
mtx(3)=readsPar([prefix,'.pro'],'lImagesPerSlab');
    nref=3;
  
    ndiff=readsPar([prefix,'.pro'],'lDiffDirections');
    nchan=readsPar([prefix,'.pro'],'lRxChannelConnected');
    nchan=length(nchan);
    loneslice=(mtx(1)*4+32)*(mtx(2)+nref)*nchan;
 

    %
     % d = fread(fid,loneslice*7*64,'*float32');
      
      %d=reshape(d,[mtx(1)*4+32,nchan,mtx(2)+nref,64,ndiff+1]);
      
    %
    Data=zeros(loneslice,(ndiff+1)*2,'single');

for i=1:(ndiff+1)*2

    fseek(fid,loneslice*(slice-1)*4,0);
    
    Data(:,i) = fread(fid,loneslice,'*float32');
    
     fseek(fid,loneslice*(mtx(3)-slice)*4,0);
     
end

    

    
 fclose(fid);
 
 
 Data=reshape(Data,[mtx(1)*4+32,nchan,mtx(2)+nref,2*(ndiff+1)]);
 
 Data=permute(Data,[1,3,4,2]);
 Data(:,3,:,:)=[];
 Data=Data(33:2:end-1,:,:,:)+1i*Data(34:2:end,:,:,:);
 
 Data(:,2:2:end,:,:)=flip(Data(:,2:2:end,:,:),1);
 
 Data=phase_corr(Data(:,3:end,:,:),Data(:,1:2,:,:));
%  img=fft2c(Data);
%  img2=sqrt(sum(abs(img).^2,4));
%  figure;imshow(abs(img2(:,:,1)),[]);
 fprintf('Read data took %s s\n',toc);

 
 function d1bk=phase_corr(d1b,ref)
    d1b=ifft1c(d1b,1);
    ref=ifft1c(ref,1);
    ph=angle(ref);
 
    d1b(:,1:2:end,:,:,:)=d1b(:,1:2:end,:,:,:).*repmat(exp(-1i*ph(:,1,:,:,:)),[1,ceil(size(d1b,2)/2),1,1,1]);
    d1b(:,2:2:end,:,:,:)=d1b(:,2:2:end,:,:,:).*repmat(exp(-1i*ph(:,2,:,:,:)),[1,floor(size(d1b,2)/2),1,1,1]);
     d1bk=fft1c(d1b,1);

     

