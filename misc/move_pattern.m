function res=move_pattern(pat,dest)

a=find_pattern(pat);

for i=1:length(a)
   movefile(a{i},dest); 
end


