function ecg_plot(fpattern)
% ecg_normalize(fpattern)

dir_str = dir(fpattern);

for i=1:length(dir_str)
    a = load(dir_str(i).name);
    

figure;
plot(a);
xlim([0,length(a)]);

set(gcf,'Units','pixels','Position',[6,653,1590,461]);

t = strrep(dir_str(i).name,'_','\_');
title(t);
end
 
