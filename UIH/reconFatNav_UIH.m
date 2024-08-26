function reconFatNav_UIH(dname)
%exclude the reference scan if data is undersampled;
[n,suf]=strtok(dname,'_');

fname=name4pat(fullfile(dname,['*',suf,'.raw']),1);
xmlname=name4pat(fullfile(dname,['*',suf,'.prot']),1);
xml=parseXML(xmlname);
if isempty(fname) disp('no file found');return; end
try
    [~,ppadata,phasecor,feed_back] = Read_UIH_Raw_v5_6(fname,1);
catch
    [~,ppadata,phasecor,feed_back] = Read_UIH_Raw_v5_6_for_fatnav(fname,1);
end
%rawdata = Read_UIH_Raw_v5_5(fname,1);
%d=ReadLineDataV1(fname,1);
%%

fb=permute(feed_back,[1,5,4,2,3]); %[ro,ch,shots,par,pe]
clear feed_back;

fb=ifft1c(fb,1);
fb=fb(end/4+1:end-end/4,:,:,:,:);
   
%%

shift=readPar_uih(xml,'FatNavFOVShift'); %shifts along [ro, pe, par]
voxsize=readPar_uih(xml,'FatNavVoxSize');
fov=readPar_uih(xml,'FatNavFOVSize'); %FOV along [ro, pe, par]

af=readPar_uih(xml,'FatNavPATFactor');

sz=size(fb);

sz(2:3)=fov(2:3)/voxsize;

d=zeros(sz);

d(:,1:size(fb,2),1:size(fb,3),:,:)=fb;

    
inc_row=1:size(fb,2);
inc_col=1:size(fb,3);
clear fb;
%%

if af==1
  
fd=ifft1c(ifft1c(d,2),3);
im=squeeze(sos(fd,4));

else 

 
 
lRefLine=readPar_uih(xml,'FatNavPATRefLines');

dacs=d(:,end/2-lRefLine/2+1:end/2+lRefLine/2,end/2-lRefLine/2+1:end/2+lRefLine/2,:,1);

im=zeros(sz(1)/2,sz(2),sz(3),sz(5)-1,'single');

dref=d(:,:,:,:,1);
smask=sos(d(:,:,:,:,2),4)>0;
save([dname,'_smask.mat'],'smask');
save([dname,'_ref.mat'],'dref');
tic;
    for j=1:size(dacs,1)
        tmp=shiftdim(d(j,:,:,:,2:end),1);    
 
     
        kres=0*tmp;
       
        kres(inc_row,inc_col,:,:)=GRAPPA_new(tmp(inc_row,inc_col,:,:),squeeze(dacs(j,:,:,:)),[7,7],0.01,[af,af]);
     
        im(j,:,:,:) = sos(ifft1c(ifft1c(kres,1),2),3); %need to check before
        time_left(j,size(dacs,1),toc);
    end

end

%%

    o.voxsize=voxsize*[1,1,1];

    pos(1)=shift(2);
    pos(2)=-shift(3);
    pos(3)=-shift(1);


    o.center=pos; %coord in dicom convention
    
    rotmat=[0,0,-1;1,0,0;0,-1,0]';
    
    o.rotmat=rotmat;
  
    o.sz=sz;
    o=center2pos_o(o,o.center);
        
    mat2niigz(im,[],[dname,'_fatnav.nii.gz'],true,o);
    

%the shift along readout seems incorrrect.