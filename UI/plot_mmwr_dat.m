
function plot_mmwr_dat(dname)

fname=name4pat(fullfile(dname,'*.mmwr.dat'),1);
fid = fopen(fname);

C = textscan(fid, '%s',18);
d = textscan(fid, '%f,%f,%f',Inf);

fclose(fid);

      %%
figure;

plot(d{1}/100000,d{2},'r-');






