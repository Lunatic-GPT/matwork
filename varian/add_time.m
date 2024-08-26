function t =add_time(t1,t2)
%t = add_time(fid_prefix)
% t1 and t2 are a 1*3 vectors containing hr, min, and sec.

tmp=t1+t2;
t=zeros(1,3);
t(3)=mod(tmp(3),60);
tmp(2)=tmp(2)+(tmp(3)-t(3))/60;

t(2)=mod(tmp(2),60);
t(1)=tmp(1)+(tmp(2)-t(2))/60;
