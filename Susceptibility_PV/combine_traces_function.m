function combine_traces_function(dname)
% input parameter
% dname: folder name where the traces of individual paths are saved
% Output: a file named dname.swc for the combined paths.

dmv=dir(fullfile(dname,'*.swc'));

fid=fopen(sprintf('%s.swc',dname),'w');
offset=0;
for j=1:length(dmv)
    pos=read_swc(fullfile(dmv(j).folder,dmv(j).name));
    
    
    write_swc(fid,pos,offset);
    
    offset=offset+size(pos,1);
end
fclose(fid);



function pos=read_swc(fname)

fid=fopen(fname,'r');

a=textscan(fid,'%f %f %f %f %f %f %f','CommentStyle','#');
fclose(fid);
pos=zeros(length(a{1}),3);

for i=1:length(a{1})
    pos(i,:)=[a{3}(i),a{4}(i),a{5}(i)];
end


function pos=write_swc(fid, pos,offset)

%: the first point number will be offset+1


for i=1:size(pos,1)
    if i==1
        fprintf(fid,'%d 0 %f %f %f 0 -1\n',i+offset,pos(i,:));
    else
        fprintf(fid,'%d 0 %f %f %f 0 %d\n',i+offset,pos(i,:),offset+i-1);
    end
end





