
function [rows,cols]=fullySampledRegion(mref)

m2=clusterize2(mref,10);
s2=sum(m2>0,2);
s1=sum(m2>0,1);

rows = s2>=max(s2)-2;
cols=s1>=max(s1)-2;
