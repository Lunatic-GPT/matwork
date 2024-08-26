function add_rdk_subject(sid)


if ~exist('sid','var')
    i = 3;
elseif strcmp(sid,'Zong')
    i= 1;    
elseif strcmp(sid,'Peng')
     i=2;
elseif strcmp(sid,'Hua')
     i=3;
elseif strcmp(sid,'K004')
     i=4;
else
    error('unknown sid');
end


root = '/home/zong/rdk';
flist = {'Zong/rdk_12_01_08/afni','Peng/RDK_12_16_2008/afni','Hua/rdk_12_04_08/afni','KELIRIS004_GeorgeSister/rdk/afni','KELIRIS003/afni','KELIRIS005/rdk/afni'};

cd(fullfile(root,flist{i}));
