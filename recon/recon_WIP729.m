%function a=recon_WIP729(fname)
fname='meas_MID89_tse_UHF_WIP729C_obl_TE28ms_FID3576.dat';
a=readMeasDat(fname);

%%
a=reshape(a(:,33:end),[1408,32,8,30,32]);
%%
fid=fopen('meas_MID89_tse_UHF_WIP729C_obl_TE28ms_FID3576.dat.log');
c=textscan(fid,'%s');
fclose(fid);
k=[];
for i=1:length(c{1})
    
    if strcmp(c{1}{i},'Slice')
        k(end+1)=str2num(c{1}{i+3});
    end
end
k=k(2:end);

k=reshape(k,[8,30,32]);
k=k(:,1,:);

%%
img=zeros(1408,704,6,32);
for sl=12:17;
b=reshape(a(:,:,:,sl,:),[1408,32,256]);

b2=permute(b,[1,3,2]);
b3=zeros(1408,704,32);
b3(:,k(:),:)=b2;

%Rk2=zeros( 1408,704,30,  32);
%Rk3=Rk2;


    
 % Rk2(:,:,i,:)=myGRAPPA_xp(b2,k(:),337:368,[-4,-1,2,5],-1:1,0);
  %Rk3(:,:,i,:)=myGRAPPA_xp(squeeze(b2(:,:,i,:)),k(:),337:368,[-5,-2,1,4],-1:1,0);
  Rk=myGRAPPA_xp(b3,k(:),337:368,[-4,-1,2,5],-1:1,0);
  ind=3:3:704;
  b3(:,ind,:)=Rk(:,ind,:);
  
  Rk=myGRAPPA_xp(b3,k(:),337:368,[-5,-2,1,4],-1:1,0);
  ind=1:3:704;
  b3(:,ind,:)=Rk(:,ind,:);
  img(:,:,sl-11,:)=fft2c(b3);
end

writeanalyze(abs(img),'recon_89_6slice',[0.3,0.3,0.8]);





