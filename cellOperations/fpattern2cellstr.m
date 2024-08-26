function name=fpattern2cellstr(fpat,scan)


for i=1:length(scan)
 name{i}=sprintf(fpat,scan(i));
end