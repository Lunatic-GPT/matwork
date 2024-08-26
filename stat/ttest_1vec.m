function [p,t] = ttest_1vec(x)
% [p,t] = ttest_1vec(x)
t = mean(x)/std(x)*sqrt(length(x));

if isnan(t)
    t = 0;
end
p = tTest(length(x)-1,t);