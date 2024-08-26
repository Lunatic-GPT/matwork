function rdk2_gen_seq(pattern)

sdir = dir(pattern);
d = fileparts(pattern);
nf = length(sdir);
s = zeros(nf,8);
for i=1:nf
    load(fullfile(d,sdir(i).name));
    cond = params.condlist;
    cond_i = cond(2:13:end);
    s(i,:) = ceil(cond_i/2)-1;
end

save_mat_int(s,'stim_seq.1D');

    
    
