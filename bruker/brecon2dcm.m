function brecon2dcm(d,irecon)

if ~exist('irecon','var')
    irecon=1;
end

im=read_2dseq(d,irecon);
nr=readbPar(fullfile(d,'acqp'),'NR',true);
im=reshape(im,[size(im,1),size(im,2),size(im,3)/nr,nr]);

data=im;
if irecon==1
save(d,'data');
else
    save([d,'_',num2str(irecon),'.dcm'],'data');
end


