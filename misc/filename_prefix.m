function res=filename_prefix(fname,prefix)

[d,f,a]=fileparts(fname);

res=fullfile(d,[prefix,f,a]);