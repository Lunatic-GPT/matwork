function res=normalizedResidualSD(d,Rsquare,npar)
% residual standard deviation after a linear fit normalized by the data mean
% dof: residual degrees of freedom, length(d)-n_fit_par-1;

dof=length(d)-npar-1;
m=mean(d);
TSS= sum((d-m).^2);
res = (1-Rsquare)*TSS/dof/m/m;
res=sqrt(res);