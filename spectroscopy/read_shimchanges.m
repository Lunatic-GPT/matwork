function change=read_shimchanges(user)

if ~exist('user','var')
    user='Xiaopeng';
end

if exist('/home/xiaopeng/labhome/vnmrsys/fastmap','dir') 
    file = '/home/xiaopeng/labhome/vnmrsys/fastmap/change';
else
    file='/data/users/xiaopeng/vnmrsys/fastmap/change';
end

if strcmp(user,'tao') && exist('/home/xiaopeng/tao/fastmap','dir')
    file = '/home/xiaopeng/tao/fastmap/change';
elseif strcmp(user,'tao') && exist('/data/users/tao/vnmrsys/fastmap','dir')
    file = '/data/users/tao/vnmrsys/fastmap/change';
end

d=textread(file,'%s');


change=zeros(1,8);

for i=1:8
    change(i)=str2double(d{i+1});
end

tmp=change(4);
change(4)=change(5);
change(5)=tmp;
%sv(1)=readPar('procpar','x1');
%sv(2)=readPar('procpar','y1');
%sv(3)=readPar('procpar','z1c');
%sv(4)=readPar('procpar','xz');
%sv(5)=readPar('procpar','yz');
%sv(6)=readPar('procpar','z2c');
%sv(7)=readPar('procpar','xy');
%sv(8)=readPar('procpar','x2y2');