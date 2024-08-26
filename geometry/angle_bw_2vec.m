function res=angle_bw_2vec(a,b)
% results in degrees
% a is [n1,n2,n3,..,3];
% b is [1,3];

nd=ndims(a);

sz1=ones(1,nd);
sz1(end)=3;

b=reshape(b,sz1);

sz2=size(a);
sz2(end)=1;
b=repmat(b,sz2);


d=sum(a.*b,nd);

res=acos(d./sos(a,nd)./sos(b,nd))*180/pi;




