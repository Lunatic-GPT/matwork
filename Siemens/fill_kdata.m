function d2=fill_kdata(data,lin,par)

nro=size(data,1);
lin=lin-min(lin);
par=par-min(par);

lin_unq=unique([lin(:),par(:)],'rows');

lin=lin+1;
par=par+1;

nlin=max(lin);
npar=max(par);

nrep=size(lin(:),1)/size(lin_unq,1);

d2=zeros(nro,nrep,nlin,npar);
irep=ones(nlin,npar);

for i=1:length(lin)
    if mod(i,1000)==0
        disp(i);
    end
    d2(:,irep(lin(i),par(i)),lin(i),par(i))=data(:,i);
    irep(lin(i),par(i))=irep(lin(i),par(i))+1;
end

d2=permute(d2,[1,3,4,2]);





