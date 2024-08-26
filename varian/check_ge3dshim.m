function check_ge3dshim()
home = '/home/xiaopeng/labhome';
shim = fullfile(home,'vnmrsys/shims');
gshimdir = fullfile(home,'vnmrsys/gshimdir/data');
cd(gshimdir);

%recon_ge3dshim('B0.1.fid');
%recon_ge3dshim('B0.A.fid');
%recon_ge3dshim('B0.B.fid');
%read_fdf('B0.fdf','B0.A+orig');
%read_fdf('B0.mask.fdf','B0.A+orig');

[err,info] = BrikInfo('B0.A+orig');

m=read_fdf('B0.mask.fdf');

coils={'x1','y1','z1c','xz','yz','z2c','xy','x2y2'};
offset = [1600,1600,1600,1920,1920,1600,1920,1920];
val=zeros(1,8);

for i=1:8
 val(i)=readPar([shim,'/start0'],coils{i});
 
end
c=read_gshimout(fullfile(gshimdir,'gshim.out'),coils);

val=(val+c)./offset;

calib=BrikLoad(fullfile(home,'vnmr/gshimdir/calib/shimmap.12cmHD.40.40.40+orig'));
sz=size(calib);
val=shiftdim(val,-2);
val=repmat(val,[sz(1:3),1]);
shimmap=val.*calib;
shimmap=sum(shimmap,4);
history = 'check_ge3dshim';

WriteBrikEZ(shimmap,info,history,'shimmap');
B0=read_fdf('B0.fdf');

figure;
x=shimmap(m>0);
y=B0(m>0);
plot(x,y,'o');
[cc,p]=corrcoef(x,y);
disp(cc(1,2));
disp(p(1,2));

function change=read_gshimout(fname,coils)

%a=textread(parfile,'%s%s%s%s%s%s%s%s%s%s%s',);
fid=fopen(fname);
a=textscan(fid,'%s');
change=zeros(1,length(coils));
for i=1:length(coils)
ind = strmatch(coils{i},a{1},'exact');
if length(ind) ~=1
    error('Field not found');
end

change(i)=str2double(a{1}{ind+1});
end






