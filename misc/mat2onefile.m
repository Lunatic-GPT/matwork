function mat2onefile(prefix,ind,dim)

for i=1:length(ind)

dir_str=dir(sprintf('%s_%d_*.mat',prefix,ind(i)));
disp(i);
disp(length(dir_str));
if i==1
im=ri(dir_str.name);
else
tmp=ri(dir_str.name);

im=cat(dim,im,tmp);

end

end
save([prefix,'_all.mat'],'im');