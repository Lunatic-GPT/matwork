function movefile_pattern(pattern,dest)
% movefile_pattern(pattern,dest)
src=find_pattern(pattern);
src=str2cell(src);
for i=1:length(src)
 movefile(src{i},dest);

end