function recon_tfl(filename,varargin)
% Example:
%
% kSpaceData = readMeasDat('meas_yyT2.dat',[204 256 70]);
% kSpaceData = readMeasDat('meas_yyT2.dat',[204 256 70],'Echo',2,'Rep',50);


[d,kindex]=readMeasDat(filename,'coilScaling.txt',[]);
%%
prefix=strtok(filename,'.');

 
mtx(1)=readsPar([prefix,'.pro'],'lBaseResolution');
mtx(2)=readsPar([prefix,'.pro'],'lPhaseEncodingLines');
mtx(3)=readsPar([prefix,'.pro'],'lConc');
%na= readsPar([prefix,'.pro'],'lAverages');   
    nchan=readsPar([prefix,'.pro'],'lRxChannelConnected');
    nchan=length(nchan);
    nrep=readsPar([prefix,'.pro'],'lRepetitions');
    nrep=nrep+1;
    
  %%  
  %  loneslice=(mtx(1)*4+32)*(mtx(2)+nref)*nchan;
 
 k=kindex(1:nchan:end);
 k=k(2:end);
 lk=length(k);
 %%
 d2=0.5*(d(1:2:end,:)+d(2:2:end,:));
 d2_noise=d2(:,1:nchan);
 d2=d2(:,nchan+1:end);
d2=reshape(d2,[mtx(1),nchan,lk/nrep,nrep,mtx(3)]);

d2=permute(d2,[1,3,5,2,4]);
%%
 k2=k(1:lk/nrep);
 
 d3=zeros([mtx(1:3),nchan,nrep]);
 d3(:,k2,:,:,:)=d2;
 %%
 [im,fd]=recon_grappa2D(d2,k2,mtx(2));
 
 
 save(['recon_',prefix],'im','fd');
    
    
 
 