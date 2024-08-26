function [p,t]=ttest_2vec(d1,d2)
%[p,t]=ttest_2vec(d1,d2)
% d1, d2's are two groups.  They can have different lengths.
m1 = mean(d1);
m2 = mean(d2);
n1 = length(d1);
n2 = length(d2);

s2 = (var(d1,1)*n1+var(d2,1)*n2)/(n1+n2-2);
t = (m1-m2)/sqrt(s2*(1/n1+1/n2));
if ~isnan(t)
p = tTest(n1+n2-2,t);
else
    p = 0;
end