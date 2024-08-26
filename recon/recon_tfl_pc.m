function recon_tfl_pc(filename,varargin)
% Example:
%
% kSpaceData = readMeasDat('meas_yyT2.dat',[204 256 70]);
% kSpaceData = readMeasDat('meas_yyT2.dat',[204 256 70],'Echo',2,'Rep',50);

if ~exist([filename,'_protocol'],'file')
  
  checkMeasDat(filename);
end


prefix=strtok(filename,'.');
fid = fopen(filename,'rb');
nProtHeaderLen = double(fread(fid,1,'*int32'));


fseek(fid, nProtHeaderLen - 4, 0);
 
mtx(1)=readsPar([prefix,'.pro'],'lBaseResolution');
mtx(2)=readsPar([prefix,'.pro'],'lPhaseEncodingLines');
mtx(3)=readsPar([prefix,'.pro'],'lConc');
na= readsPar([prefix,'.pro'],'lAverages');   
    nchan=readsPar([prefix,'.pro'],'lRxChannelConnected');
    nchan=length(nchan);
  %  loneslice=(mtx(1)*4+32)*(mtx(2)+nref)*nchan;
 
data=fread(fid,'*float');
%data=fread(fid,'*int64');

%data=double(data);
    k=load(sprintf('%s.dat_kindex.log',prefix));
    k(1,:)=[];
    lk=size(k,1);
 fclose(fid);
 
 
 off=(mtx(1)*4+32)*nchan;
m2=[mtx(1)*4+32,nchan,2*na,lk/2/na,mtx(3)];

 data=data(off:prod(m2)+off-1);
 
 data=reshape(data,m2);
 data=data(33:end,:,:,:,:);
 
 data=permute(data,[1,4,5,2,3]);
  data=data(1:2:end,:,:,:,:)+1i*data(2:2:end,:,:,:,:);
 data(1,:,:,:,:)=0;
 
 k2=k(1:2*na:end,2);

   %[im,fd]=recon_grappa2D(data,k2,mtx(2));
 
 [im,fd]=recon_grappa2D_Lustig(data,k2,mtx(2));
 
 save(prefix,'im','fd');
    
    
 
 