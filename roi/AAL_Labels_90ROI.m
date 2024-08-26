function region=AAL_Labels_90ROI

fid=fopen('atlas90-label-description.txt','r');
res=textscan(fid,'%d%s%d');
region=cell(1,45);

for i=1:45
   region{i}=res{2}{2*i}(1:end-2);  
end
fclose(fid);