function reconTSE_UIH(dname)

[n,suf]=strtok(dname,'_');

fname=name4pat(fullfile(dname,['*',suf,'.raw']),1);
xmlname=name4pat(fullfile(dname,['*',suf,'.prot']),1);
xml=parseXML(xmlname);
if isempty(fname) disp('no file found');return; end
[data,ppadata,phasecor,feed_back] = Read_UIH_Raw_v5_6(fname);

data=squeeze(data); %[pe,par,ro,ch]
data=permute(data,[3,1,2,4]);%change to [ro,pe,par,ch]

fd=ifft1c(data,1);
fd=fd(end/4+1:end-end/4,:,:,:);
   
sz=size(fd);
img=zeros(sz(1:3));
for i=1:size(fd,1)

    

    img(i,:,:)=tmp;

end

