function fdf2afni(d)


d2=[d,'.fid'];
d1=[d,'.img'];
ns = readPar(d2,'ns');
nread=readPar(d2,'nread');
nphase=readPar(d2,'nphase');
nread=nread/2;

str=dir(fullfile(d1,'*.fdf'));
if mod(length(str),ns)~=0
    error('number of files wrong');
end

ad=length(str)/ns;

b=zeros(nread,nphase,ns,ad);
for i=1:ns
    for j=1:ad
    
     name = sprintf('slice%03dimage%03decho001.fdf',i,j);

     tmp=read_fdf(fullfile(d1,name));
     b(:,:,i,j)=squeeze(tmp);
    end
end

lro=readPar(d2,'lro');
lpe=readPar(d2,'lpe');

pro=readPar(d2,'pro');
ppe=readPar(d2,'ppe');
pss0=readPar(d2,'pss0');

thk=readPar(d2,'thk');
orient=readPar(d2,'orient');
orient=orient(2:end-1);

lro=lro*10;
lpe=lpe*10;

orig_delta=[pro,ppe,pss0;lro/nread,lpe/nphase,thk];
[b,orig_delta]=reorient_data(b,orig_delta,orient);

info.ORIGIN = orig_delta(1,:);
info.DELTA = orig_delta(2,:);

write_afni(b,info,d);


