function [d,ind_rm]=unique_1col(d,icol)

dd=diff(d(:,icol));

ind=find(dd==0);

d(ind+1,:)=[];

ind_rm=ind+1;

[tmp,ind]=sort(d(:,icol));
d=d
