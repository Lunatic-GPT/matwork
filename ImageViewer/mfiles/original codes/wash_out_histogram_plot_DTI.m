mu = -0.1459; % mean of malignant distribution

ADC_vec_b= [];
ADC_vec_m = [];

x_hist = 0:0.05:3.5;

%close all
disp(' ---------- NEW DATA --------- ')

sourcy = 'B298964';
source = ['/export/res/breast/' sourcy '/research_results_DTI_analysis_model_set_up/slice3/ADC_vector.mat'];
load(source, 'ADC_lesion_vec_shifted');
ADC_vec_b = [ADC_vec_b ADC_lesion_vec_shifted];
temp = ADC_lesion_vec_shifted(ADC_lesion_vec_shifted~=0);
meany = mean(temp)-1

sourcy = 'B162582';
source = ['/export/res/breast/' sourcy '/research_results_DTI_analysis_model_set_up/slice3/ADC_vector.mat'];
load(source, 'ADC_lesion_vec_shifted');
ADC_vec_b = [ADC_vec_b ADC_lesion_vec_shifted];
temp = ADC_lesion_vec_shifted(ADC_lesion_vec_shifted~=0);
meany = mean(temp)-1

sourcy = 'B153241';
source = ['/export/res/breast/' sourcy '/research_results_DTI_analysis_model_set_up/slice2/ADC_vector.mat'];
load(source, 'ADC_lesion_vec_shifted');
ADC_vec_b = [ADC_vec_b ADC_lesion_vec_shifted];
temp = ADC_lesion_vec_shifted(ADC_lesion_vec_shifted~=0);
meany = mean(temp)-1

sourcy = 'B49217';
source = ['/export/res/breast/' sourcy '/research_results_DTI_analysis_model_set_up/slice2/ADC_vector.mat'];
load(source, 'ADC_lesion_vec_shifted');
ADC_vec_b = [ADC_vec_b ADC_lesion_vec_shifted];
temp = ADC_lesion_vec_shifted(ADC_lesion_vec_shifted~=0);
meany = mean(temp)-1

sourcy = 'B299986';
source = ['/export/res/breast/' sourcy '/research_results_DTI_analysis_model_set_up/slice2/ADC_vector.mat'];
load(source, 'ADC_lesion_vec_shifted');
ADC_vec_b = [ADC_vec_b ADC_lesion_vec_shifted];
temp = ADC_lesion_vec_shifted(ADC_lesion_vec_shifted~=0);
meany = mean(temp)-1

sourcy = 'N309953benign';
source = ['/export/res/breast/' sourcy '/research_results_DTI_analysis_model_set_up/slice3/ADC_vector.mat'];
load(source, 'ADC_lesion_vec_shifted');
ADC_vec_b = [ADC_vec_b ADC_lesion_vec_shifted];
temp = ADC_lesion_vec_shifted(ADC_lesion_vec_shifted~=0);
meany = mean(temp)-1


sourcy = 'M298821';
source = ['/export/res/breast/' sourcy '/research_results_DTI_analysis_model_set_up/slice3/ADC_vector.mat'];
load(source, 'ADC_lesion_vec_shifted');
ADC_vec_m = [ADC_vec_m ADC_lesion_vec_shifted];

sourcy = 'M301431';
source = ['/export/res/breast/' sourcy '/research_results_DTI_analysis_model_set_up/slice3/ADC_vector.mat'];
load(source, 'ADC_lesion_vec_shifted');
ADC_vec_m = [ADC_vec_m ADC_lesion_vec_shifted];

sourcy = 'M179579';
source = ['/export/res/breast/' sourcy '/research_results_DTI_analysis_model_set_up/slice3/ADC_vector.mat'];
load(source, 'ADC_lesion_vec_shifted');
ADC_vec_m = [ADC_vec_m ADC_lesion_vec_shifted];

sourcy = 'M7803';
source = ['/export/res/breast/' sourcy '/research_results_DTI_analysis_model_set_up/slice3/ADC_vector.mat'];
load(source, 'ADC_lesion_vec_shifted');
ADC_vec_m = [ADC_vec_m ADC_lesion_vec_shifted];

sourcy = 'M273920';
source = ['/export/res/breast/' sourcy '/research_results_DTI_analysis_model_set_up/slice3/ADC_vector.mat'];
load(source, 'ADC_lesion_vec_shifted');
ADC_vec_m = [ADC_vec_m ADC_lesion_vec_shifted];

sourcy = 'N84432';
source = ['/export/res/breast/' sourcy '/research_results_DTI_analysis_model_set_up/slice3/ADC_vector.mat'];
load(source, 'ADC_lesion_vec_shifted');
ADC_vec_m = [ADC_vec_m ADC_lesion_vec_shifted];

sourcy = 'M297672';
source = ['/export/res/breast/' sourcy '/research_results_DTI_analysis_model_set_up/slice3/ADC_vector.mat'];
load(source, 'ADC_lesion_vec_shifted');
ADC_vec_m = [ADC_vec_m ADC_lesion_vec_shifted];

sourcy = 'B49217b';
source = ['/export/res/breast/' sourcy '/research_results_DTI_analysis_model_set_up/slice3/ADC_vector.mat'];
load(source, 'ADC_lesion_vec_shifted');
ADC_vec_m = [ADC_vec_m ADC_lesion_vec_shifted];

sourcy = 'M303424';
source = ['/export/res/breast/' sourcy '/research_results_DTI_analysis_model_set_up/slice3/ADC_vector.mat'];
load(source, 'ADC_lesion_vec_shifted');
ADC_vec_m = [ADC_vec_m ADC_lesion_vec_shifted];

% sourcy = 'N309953';
% source = ['/export/res/breast/' sourcy '/research_results_DTI_analysis_model_set_up/slice11/ADC_vector.mat'];
% load(source, 'ADC_lesion_vec_shifted');
% ADC_vec_m = [ADC_vec_m ADC_lesion_vec_shifted];

temp = ADC_vec_m(ADC_vec_m~=0);
mal = temp;
mean_malignant = mean(temp)
std_malignant = std(temp)
figure; 
subplot(2,2,1);
%hist(temp,x_hist); title('malignant');
h_m = hist(temp,x_hist);
hn_m = h_m/(length(temp)*0.05); % bin width
bar(x_hist,hn_m);
% figure; histfit(temp,50); title('histogram incl. fit');

%ADC_mean = mean(ADC_vec_m(ADC_vec_m~=0))
%ADC_std = std(ADC_vec_m(ADC_vec_m~=0))

x = x_hist;


% d = (x-ADC_mean)/ADC_std;
% gauss = 18474*(1/sqrt(2*pi))*exp(-d.^2/2);
% gauss = (1/sqrt(2*pi))*exp(-d.^2/2)/ADC_std;
% figure; plot(x,gauss);
% hold on
% 
% h_m = hist(temp,x);
% hn_m = h_m/(length(temp)*3); % 3 is the bin width
% bar(x,hn_m);
% hold off;
% 
% figure
temp = ADC_vec_b(ADC_vec_b~=0);
% ben = temp;
mean_benign = mean(temp)
std_benign = std(temp)
% h_b = hist(temp,x);
% hn_b = h_b/(length(temp)*3); 
% bar(x,hn_b);

subplot(2,2,3);
h_m = hist(temp,x_hist);
hn_m = h_m/(length(temp)*0.05); % bin width
bar(x_hist,hn_m);
%hist(temp,x_hist); title('benign');
% figure; hist(ADC_vec_b(ADC_vec_b~=0),50); title('benign');

% figure; hist(temp,x); title('malignant distribution');
% plot(x,gauss);
% hold off



% larger=numel(ADC_vec(ADC_vec>0))
% smaller=numel(ADC_vec(ADC_vec<0))

% ADC_vec = [ADC_vec_b ADC_vec_m];
% figure; hist(ADC_vec(ADC_vec~=0),50); title('all lesions');

ADC_vec_b= [];
ADC_vec_m = [];

%close all

sourcy = 'B298964';
source = ['/export/res/breast/' sourcy '/research_results_DTI_analysis_model_set_up/slice3/ADC_vector.mat'];
load(source, 'ADC_lesion_vec_shifted');
ADC_vec_b = [ADC_vec_b ADC_lesion_vec_shifted/0.958];

sourcy = 'B162582';
source = ['/export/res/breast/' sourcy '/research_results_DTI_analysis_model_set_up/slice3/ADC_vector.mat'];
load(source, 'ADC_lesion_vec_shifted');
ADC_vec_b = [ADC_vec_b ADC_lesion_vec_shifted/1.010];

sourcy = 'B153241';
source = ['/export/res/breast/' sourcy '/research_results_DTI_analysis_model_set_up/slice2/ADC_vector.mat'];
load(source, 'ADC_lesion_vec_shifted');
ADC_vec_b = [ADC_vec_b ADC_lesion_vec_shifted/1.375];

sourcy = 'B49217';
source = ['/export/res/breast/' sourcy '/research_results_DTI_analysis_model_set_up/slice2/ADC_vector.mat'];
load(source, 'ADC_lesion_vec_shifted');
ADC_vec_b = [ADC_vec_b ADC_lesion_vec_shifted/1.439];

sourcy = 'B299986';
source = ['/export/res/breast/' sourcy '/research_results_DTI_analysis_model_set_up/slice2/ADC_vector.mat'];
load(source, 'ADC_lesion_vec_shifted');
ADC_vec_b = [ADC_vec_b ADC_lesion_vec_shifted/1.093];

sourcy = 'N309953benign';
source = ['/export/res/breast/' sourcy '/research_results_DTI_analysis_model_set_up/slice3/ADC_vector.mat'];
load(source, 'ADC_lesion_vec_shifted');
ADC_vec_b = [ADC_vec_b ADC_lesion_vec_shifted/1.031];


sourcy = 'M298821';
source = ['/export/res/breast/' sourcy '/research_results_DTI_analysis_model_set_up/slice3/ADC_vector.mat'];
load(source, 'ADC_lesion_vec_shifted');
ADC_vec_m = [ADC_vec_m ADC_lesion_vec_shifted/1.286];

sourcy = 'M301431';
source = ['/export/res/breast/' sourcy '/research_results_DTI_analysis_model_set_up/slice3/ADC_vector.mat'];
load(source, 'ADC_lesion_vec_shifted');
ADC_vec_m = [ADC_vec_m ADC_lesion_vec_shifted/1.689];

sourcy = 'M179579';
source = ['/export/res/breast/' sourcy '/research_results_DTI_analysis_model_set_up/slice3/ADC_vector.mat'];
load(source, 'ADC_lesion_vec_shifted');
ADC_vec_m = [ADC_vec_m ADC_lesion_vec_shifted/1.211];

sourcy = 'M7803';
source = ['/export/res/breast/' sourcy '/research_results_DTI_analysis_model_set_up/slice3/ADC_vector.mat'];
load(source, 'ADC_lesion_vec_shifted');
ADC_vec_m = [ADC_vec_m ADC_lesion_vec_shifted/1.021];

sourcy = 'M273920';
source = ['/export/res/breast/' sourcy '/research_results_DTI_analysis_model_set_up/slice3/ADC_vector.mat'];
load(source, 'ADC_lesion_vec_shifted');
ADC_vec_m = [ADC_vec_m ADC_lesion_vec_shifted/1.534];

sourcy = 'N84432';
source = ['/export/res/breast/' sourcy '/research_results_DTI_analysis_model_set_up/slice3/ADC_vector.mat'];
load(source, 'ADC_lesion_vec_shifted');
ADC_vec_m = [ADC_vec_m ADC_lesion_vec_shifted/1.119];

sourcy = 'M297672';
source = ['/export/res/breast/' sourcy '/research_results_DTI_analysis_model_set_up/slice3/ADC_vector.mat'];
load(source, 'ADC_lesion_vec_shifted');
ADC_vec_m = [ADC_vec_m ADC_lesion_vec_shifted/1.378];

sourcy = 'B49217b';
source = ['/export/res/breast/' sourcy '/research_results_DTI_analysis_model_set_up/slice3/ADC_vector.mat'];
load(source, 'ADC_lesion_vec_shifted');
ADC_vec_m = [ADC_vec_m ADC_lesion_vec_shifted/1.385];

sourcy = 'M303424';
source = ['/export/res/breast/' sourcy '/research_results_DTI_analysis_model_set_up/slice3/ADC_vector.mat'];
load(source, 'ADC_lesion_vec_shifted');
ADC_vec_m = [ADC_vec_m ADC_lesion_vec_shifted/1.301];

% sourcy = 'N309953';
% source = ['/export/res/breast/' sourcy '/research_results_DTI_analysis_model_set_up/slice11/ADC_vector.mat'];
% load(source, 'ADC_lesion_vec_shifted');
% ADC_vec_m = [ADC_vec_m ADC_lesion_vec_shifted];

temp = ADC_vec_m(ADC_vec_m~=0);
mal = temp;
mean_malignant_rel = mean(temp)-1
std_malignant_rel = std(temp)
subplot(2,2,2); hist(temp-1,x_hist-1); title('malignant relative');
axis([-1.5 1.5 0 2])
% figure; histfit(temp,50); title('histogram incl. fit');

ADC_mean = mean(ADC_vec_m(ADC_vec_m~=0))-1;
ADC_std = std(ADC_vec_m(ADC_vec_m~=0));

x = x_hist-1;

d = (x-ADC_mean)/ADC_std;

% gauss = (1/sqrt(2*pi))*exp(-d.^2/2)/ADC_std;
% subplot(2,2,2); plot(x,gauss);
% hold on

h_m = hist(temp-1,x);
hn_m = h_m/(length(temp)*0.05); % bin width
bar(x,hn_m);
hold on;
gauss = (1/sqrt(2*pi))*exp(-d.^2/2)/ADC_std;
subplot(2,2,2); plot(x,gauss);

hold off;

% %figure
temp = ADC_vec_b(ADC_vec_b~=0);
% ben = temp;
mean_benign_rel = mean(temp)-1
std_benign_rel = std(temp)
% h_b = hist(temp,x);
% hn_b = h_b/(length(temp)*3); 
% bar(x,hn_b);

subplot(2,2,4); 
h_m = hist(temp-1,x);
hn_m = h_m/(length(temp)*0.05); % bin width
bar(x,hn_m);
%hist(temp-1,x_hist-1); title('benign relative');
axis([-1.5 1.5 0 2.1])
% figure; hist(ADC_vec_b(ADC_vec_b~=0),50); title('benign');

% figure; hist(temp,x); title('malignant distribution');
% plot(x,gauss);
% hold off



% larger=numel(ADC_vec(ADC_vec>0))
% smaller=numel(ADC_vec(ADC_vec<0))

% ADC_vec = [ADC_vec_b ADC_vec_m];
% figure; hist(ADC_vec(ADC_vec~=0),50); title('all lesions');