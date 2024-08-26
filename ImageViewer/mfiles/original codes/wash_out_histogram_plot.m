wash_out_vec_b= [];
wash_out_vec_m = [];

sourcy = 'B298964';
source = ['/export/res/breast/' sourcy '/research_results/slice3/wash_out_vector.mat'];
load(source, 'wash_out_ROI_vec');
wash_out_vec_b = [wash_out_vec_b wash_out_ROI_vec];

sourcy = 'B162582';
source = ['/export/res/breast/' sourcy '/research_results/slice3/wash_out_vector.mat'];
load(source, 'wash_out_ROI_vec');
wash_out_vec_b = [wash_out_vec_b wash_out_ROI_vec];

sourcy = 'B153241';
source = ['/export/res/breast/' sourcy '/research_results/slice2/wash_out_vector.mat'];
load(source, 'wash_out_ROI_vec');
wash_out_vec_b = [wash_out_vec_b wash_out_ROI_vec];

sourcy = 'B49217';
source = ['/export/res/breast/' sourcy '/research_results/slice2/wash_out_vector.mat'];
load(source, 'wash_out_ROI_vec');
wash_out_vec_b = [wash_out_vec_b wash_out_ROI_vec];

sourcy = 'B299986';
source = ['/export/res/breast/' sourcy '/research_results/slice2/wash_out_vector.mat'];
load(source, 'wash_out_ROI_vec');
wash_out_vec_b = [wash_out_vec_b wash_out_ROI_vec];


sourcy = 'M298821';
source = ['/export/res/breast/' sourcy '/research_results/slice3/wash_out_vector.mat'];
load(source, 'wash_out_ROI_vec');
wash_out_vec_m = [wash_out_vec_m wash_out_ROI_vec];

sourcy = 'M301431';
source = ['/export/res/breast/' sourcy '/research_results/slice3/wash_out_vector.mat'];
load(source, 'wash_out_ROI_vec');
wash_out_vec_m = [wash_out_vec_m wash_out_ROI_vec];

sourcy = 'M179579';
source = ['/export/res/breast/' sourcy '/research_results/slice3/wash_out_vector.mat'];
load(source, 'wash_out_ROI_vec');
wash_out_vec_m = [wash_out_vec_m wash_out_ROI_vec];

sourcy = 'M7803';
source = ['/export/res/breast/' sourcy '/research_results/slice3/wash_out_vector.mat'];
load(source, 'wash_out_ROI_vec');
wash_out_vec_m = [wash_out_vec_m wash_out_ROI_vec];

sourcy = 'M273920';
source = ['/export/res/breast/' sourcy '/research_results/slice3/wash_out_vector.mat'];
load(source, 'wash_out_ROI_vec');
wash_out_vec_m = [wash_out_vec_m wash_out_ROI_vec];

sourcy = 'N84432';
source = ['/export/res/breast/' sourcy '/research_results/slice3/wash_out_vector.mat'];
load(source, 'wash_out_ROI_vec');
wash_out_vec_m = [wash_out_vec_m wash_out_ROI_vec];

sourcy = 'M297672';
source = ['/export/res/breast/' sourcy '/research_results/slice3/wash_out_vector.mat'];
load(source, 'wash_out_ROI_vec');
wash_out_vec_m = [wash_out_vec_m wash_out_ROI_vec];

sourcy = 'B49217b';
source = ['/export/res/breast/' sourcy '/research_results/slice3/wash_out_vector.mat'];
load(source, 'wash_out_ROI_vec');
wash_out_vec_m = [wash_out_vec_m wash_out_ROI_vec];

sourcy = 'M303424';
source = ['/export/res/breast/' sourcy '/research_results/slice3/wash_out_vector.mat'];
load(source, 'wash_out_ROI_vec');
wash_out_vec_m = [wash_out_vec_m wash_out_ROI_vec];

temp = wash_out_vec_m(wash_out_vec_m~=0);
% figure; hist(temp,50); title('malignant');
% figure; histfit(temp,50); title('histogram incl. fit');

wash_out_mean = mean(wash_out_vec_m(wash_out_vec_m~=0));
wash_out_std = std(wash_out_vec_m(wash_out_vec_m~=0));

x = -90:3:90;

d = (x-wash_out_mean)/wash_out_std;
% gauss = 18474*(1/sqrt(2*pi))*exp(-d.^2/2);
gauss = (1/sqrt(2*pi))*exp(-d.^2/2)/wash_out_std;
figure; plot(x,gauss);

hold on;
h_m = hist(temp,x);
hn_m = h_m/(length(temp)*3); % 3 is the bin width
bar(x,hn_m);

figure
temp = wash_out_vec_b(wash_out_vec_b~=0);
h_b = hist(temp,x);
hn_b = h_b/(length(temp)*3); 
bar(x,hn_b);
%figure; hist(wash_out_vec_b(wash_out_vec_b~=0),50); title('benign');

% figure; hist(temp,x); title('malignant distribution');
plot(x,gauss);
hold off



% larger=numel(wash_out_vec(wash_out_vec>0))
% smaller=numel(wash_out_vec(wash_out_vec<0))

% wash_out_vec = [wash_out_vec_b wash_out_vec_m];
% figure; hist(wash_out_vec(wash_out_vec~=0),50); title('all lesions');
