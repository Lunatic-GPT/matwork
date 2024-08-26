function run_dcm2dir_useHeader
a=dir2('*');
for i=1:length(a)
    
    cd(a(i).name);
    dcm2dir_useHeader(fullfile(a(i).folder,a(i).name),1);
    cd('..');
end
