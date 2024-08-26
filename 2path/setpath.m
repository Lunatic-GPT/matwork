dname={'CoordinateTransform','ImageViewer','misc','../matlab/nifti_tools','roi','stat','mat_orientation','cellOperations',...
    'Siemens','plotTools','recon','network','PVS','FatNav','.','../matlab/afni_matlab','scanhelper'};

cur_dir=pwd;
if ispc
 root=tood;
else
 root='~';   
end

for i=1:length(dname)
addpath(genpath(fullfile(root,'matwork',dname{i})));
end


dname={'espirit/espirit/utils','stat'};
for i=1:length(dname)
disp(dname{i});
addpath(genpath(fullfile(root,'matlab',dname{i})));
end

cd(root);
cd matlab/fessler
setup

cd(cur_dir);