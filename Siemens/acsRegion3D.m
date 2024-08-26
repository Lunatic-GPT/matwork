
function [acs_row,acs_col]=acsRegion3D(mref)



m2=clusterize2(mref,10);
s2=sum(m2,2);
s1=sum(m2,1);

acs_row = s2>=max(s2)-2;
acs_col=s1>=max(s1)-2;