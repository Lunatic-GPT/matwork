function ts_mean_afni(fpattern,out_prefix)
% ts_mean_afni(fpattern,out_prefix)

a = dir(fpattern);
files = '';
for i=1:length(a)
    files=sprintf('%s %s',files,a(i).name);
end
    cmd = sprintf('3dMean -prefix %s %s',out_prefix,files);
    unix(cmd);