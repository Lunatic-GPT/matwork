function res=copy_pattern(pat,dest)

a=find_pattern(pat);

for i=1:length(a)
   copyfile(a{i},dest); 
end


