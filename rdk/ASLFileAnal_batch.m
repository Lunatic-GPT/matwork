root = 'e:\rdk2';
d = {'S2','S4','S5','S10','S13','S15'};

subd = {'rdk','ret'};

for i=1:6
    for j=1:2
        cd(fullfile(root,d{i},subd{j}));
        pwd
        dir_str = dir('*.eyd');
        if length(dir_str) == 1
            ASLfileAnal(dir_str.name);
        end
    end
end



root1 = '/media/SDMINI/rdk2';
root2 = '/home/zong/rdk2';
d = {'S2','S4','S5','S10','S13','S15'};

subd = {'rdk','ret'};
for i=1:6
    for j=1:2
        cd(fullfile(root1,d{i},subd{j}));
        pwd
        dir_str = dir('*.eyd');
       
        if length(dir_str) == 1
            
            copyfile('*.mat',fullfile(root2,d{i},subd{j}));
        end
    end
end