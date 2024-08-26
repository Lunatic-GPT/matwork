function res=find_pattern(pat)


a=dir2(pat);
res={};
for i=1:length(a)
 res{i}=fullfile(a(i).folder,a(i).name);
end


a2=dir2('*');

for i=1:length(a2)
   if a2(i).isdir
       cd(a2(i).name);
        res=[res,find_pattern(pat)];
       cd ..;
   end
    
end

if nargout==0
    for i=1:length(res)
        disp(res{i});
    end
end




