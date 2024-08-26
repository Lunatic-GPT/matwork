function r=readRFfile(name,Line2Skip)

if ~exist('Line2Skip','var')
    Line2Skip=9;
end
fid = fopen(name,'r');

for i=1:Line2Skip
tmp=fgetl(fid);
end

mag=[];
ang=[];
a=fscanf(fid,'%f %f ; (%d)');

mag=a(1:3:end);
ang=a(2:3:end);



r=mag.*exp(1i*ang/180*pi);


fclose(fid);