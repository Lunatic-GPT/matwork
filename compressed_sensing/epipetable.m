function epipetable

fid=fopen('blipscale','w');

b=load('epitable.m');

db=diff(b);
for i=1:length(db)
    fprintf(fid,'%d\n',db(i));
end
fclose(fid);