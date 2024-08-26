function brecon2analyze(d,irecon)

if ~exist('irecon','var')
    irecon=1;
end

im=read_2dseq(d,irecon);
nr=readbPar(fullfile(d,'acqp'),'NR',true);
ni=readbPar(fullfile(d,'acqp'),'NI',true);
ns=readbPar(fullfile(d,'acqp'),'NSLICES',true);
%im=reshape(im,[size(im,1),size(im,2),size(im,3)*size(im,4)/nr/ni*ns,nr*ni/ns]);
im=reshape(im,[size(im,1),size(im,2),ni/ns,ns,size(im,3)*size(im,4)/ni]);

im=permute(im,[1,2,4,3,5]);
sz=size(im);
sz(end+1:5)=1;

im=reshape(im,[sz(1:2),sz(3),sz(4)*sz(5)]);

%{
or=readbPar(fullfile(d,'method'),'PVM_ObjOrderList',true);


if length(or)==size(im,3)   
    im2=zeros(size(im));
    im2(:,:,or+1,:)=im;
else
    im2=im;
end
%}
if irecon==1
write_afni(im,d);
else
    write_afni(im,[d,'_',num2str(irecon)]);
end


