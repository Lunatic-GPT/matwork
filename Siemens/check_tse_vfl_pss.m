function check_tse_vfl_pss(fname,pc_only)


if ~exist('pc_only')
    pc_only=false;
end

if pc_only

[d,lin,par]=readMeasDat(fname,20000,0,false);

else
[d,lin,par]=readMeasDat(fname,inf,0,true);
end    
prefix=strtok(fname,'.');
nch=nChan_SiemensProt([prefix,'.pro']);
d=reshape(d,[size(d,1),nch,length(lin)/nch]);

nro=size(d,1);
lin2=lin(1:nch:end);
par2=par(1:nch:end);

nref=0;
if lin2(1)==0
  lin2(1)=[];  
 end
    
for i=1:length(lin2)  
    if lin2(i)==lin2(1)
        nref=nref+1;
    else
        break;
    end
end

if lin(1)==0
    nref=nref+1;
end
%%
figure;
res=sos(d(:,:,1:nref),2);
res=mean(squeeze(res),1);
if lin(1)==0
    plot(res(2:end));
else
    plot(res);
end

[tmp,ich]=max(sum(sum(abs(d(:,:,1:nref)),3),1));

figure;
im=squeeze(abs(d(:,ich,1:nref)));
imshow(im,[0,max(im(:))*1]);
if pc_only
    return;
end


%%


lin2=lin(1:nch:end);
par2=par(1:nch:end);
lin2=lin2(nref+1:end);
par2=par2(nref+1:end);



d2=zeros(max(lin2)-min(lin2)+1,max(par2)-min(par2)+1);
d3=zeros(nro,max(lin2)-min(lin2)+1,max(par2)-min(par2)+1);
fd=fft1c(d,1);
for i=1:length(lin2(:))
    
    d2(lin2(i)+1-min(lin2),par2(i)+1-min(par2))=squeeze(mean(sos(d(:,:,i+nref),2),1));
    d3(:,lin2(i)+1-min(lin2),par2(i)+1-min(par2))=squeeze(sos(fd(:,:,i+nref),2));
    
end

prefix=strtok(fname,'.');
save([prefix,'.mat'],'d3');
%d3=d2(d2>0);
%d3=reshape(d3,[length(unique(lin2)),length(unique(par2))]);

figure;imshow(flipud(d2),[0,0.03*max(d2(:))]);

