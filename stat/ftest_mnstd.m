function [p,f]=ftest_mnstd(xmean,xstd,n)
%[p,f]=ftest_mnstd(xmean,xstd,n)
% calculate the p value given the mean, standard deviation and number of
% samples of each group.
% NOTE: standard deviation is that of the sample, not of the sample mean.
sz1=n;
sz2=length(xmean);

   f = sz1*std(xmean)^2/mean(xstd.^2);
   
   dof1 = sz2-1;
   dof2 = sz2*(sz1-1);
   p = FTest(dof1,dof2,f);