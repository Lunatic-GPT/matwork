function par=cdp(scand)
warning off;
na=readdPar(scand, 'NumberOfAverages');
tr=readdPar(scand, 'RepetitionTime');
te=readdPar(scand, 'EchoTime');
fa=readdPar(scand, 'FlipAngle');
try
ta= readdPar(scand,'Private_0051_100a');
catch
    ta=0;
end
ps=readdPar(scand,'PixelSpacing');
try
age=readdPar(scand,'PatientAge');
catch
    age='000';
end
orient=readdPar(scand,'ImageOrientationPatient');
bd=readdPar(scand,'PatientBirthDate');
sex=readdPar(scand,'PatientSex');
try
weight=readdPar(scand,'PatientWeight');
catch
    weight=0;
end

d=dir2(scand);
extp(scand);
for i=1:length(d)
    try
        in=dicomread([scand,'\',d(i).name]);
        fprintf('Matrix Size = %d*%d\n',size(in,1),size(in,2));
        break;
    catch
        continue;
    end
    
end

try
venc=readdPar(scand,'Private_0051_1014');
fprintf('VENC = %s\n',venc);
VENC = venc;
catch
    VENC = NaN;
end
seq=readsPar([scand,'.pro'],'tSequenceFileName');
if iscell(seq)
fprintf('Sequence %s\n',seq{1});

else
    
fprintf('Sequence %s\n',seq);

end
bipolar=readsPar([scand,'.pro'],'alFree[21]');
if isempty(bipolar)
    fprintf('Bipolar Off\n');
else
    fprintf('Bipolar On\n');  
end

inv=readsPar([scand,'.pro'],'alFree[22]');  %only true for fl_fq_retroz_mb and fl_fq_retroz after 1/5/2017
if isempty(inv)
    fprintf('Inv Off\n');
else
    fprintf('Inv On\n');  
end

try 
ps(3)=readdPar(scand,'SpacingBetweenSlices');

catch
    ps(3)=0;
end
try 
    fc=readdPar(scand,'Private_0019_1011');
catch
    fc='No';
end

thk=readdPar(scand,'SliceThickness');

fprintf('NumberOfAverage = %d\n',na);
fprintf('RepetitionTime = %d ms\n',tr);
fprintf('EchoTime = %6.2f ms\n',te);
fprintf('FlipAngle = %d\n',fa);
fprintf('%s\n',ta);


fprintf('Voxel Size = %5.4f*%5.4f*%5.4f; thickness=%3.2f\n',ps(1),ps(2),ps(3),thk);
fprintf('Age = %s; Birthday = %s; Sex = %s; Weight = %4.1f kg\n',age,bd,sex,weight);
fprintf('ImageOrientationPatient= ');disp(orient');

fprintf('Flow Comp = %s\n',fc);

age=str2num(age(1:3));
sex=sex=='M'; % 1: male; 0: female
par.age=age;
par.sex=sex;
par.NA=na;
par.TR=tr;
par.TE=te;
par.FA=fa;
par.TA=ta;
par.FlowComp=fc;
par.orient=orient;
par.weight=weight;
par.VoxelSize=ps;
par.Thickness=thk;
par.Birthday=bd;
par.VENC=VENC;

name=filename_prefix(scand,'DCMPar_');
save([name,'.mat'],'par');



