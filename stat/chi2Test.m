function p = chi2Test(chi2,dof)
% p = chi2Test(chi2,dof)
cmd = sprintf('cdf -t2p fict %f %d',chi2,dof);
[err,stat] = unix(cmd);
p = str2double(stat(5:end));




  