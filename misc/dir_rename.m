function dir_rename()

a =dir('s*');
ns = length(a);


for i=1:ns
   if i==1
       name = 'loc';
   elseif i==ns
       name = 'anat';
   else
       name = sprintf('s%d',i-1);
   end
  disp(['rename ', a(i).name,' to ',name]); 
  movefile(a(i).name,name);
end
    