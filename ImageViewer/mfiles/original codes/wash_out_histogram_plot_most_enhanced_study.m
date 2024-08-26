wash_out_vec_b= [];
wash_in_vec_b = [];
wash_out_vec_m = [];
wash_in_vec_m = [];

close all

sourcy = 'B298964';
source = ['/export/res/breast/' sourcy '/research_results_combined_kinetics/slice9/wash_out_vector.mat'];
load(source, 'wash_out_ROI_vec');
wo = wash_out_ROI_vec;
wash_out_vec_b = [wash_out_vec_b wash_out_ROI_vec];

sourcy = 'B298964';
source = ['/export/res/breast/' sourcy '/research_results_combined_kinetics/slice9/wash_in_vector_uptake.mat'];
load(source, 'wash_in_ROI_vec');
wi = wash_in_ROI_vec;
wash_in_vec_b = [wash_in_vec_b wi];

wob = zeros(11,1);
wo_mean = [];
wo_mean2 = [];
wi_vec = [];
for i=60:3:90 % walk through the bins
    ind = find(wi>i & wi<i+3);
    if numel(ind) ~= 0
    for j = 1:size(ind,1)
        wo_mean2 = [wo_mean2 wo(ind(j))];
    end  
    wi_vec(((i-60)/3)+1) = i+3;
    wo_mean(((i-60)/3)+1) = mean(wo_mean2);
    wob(((i-60)/3)+1,1:size(wo_mean2,1)) = wo_mean2; 
    wo_mean2 = [];
    end
end
% figure
% plot(wi_vec,wo_mean,'d'); xlabel('wi (degree)'), ylabel('wo (degree)'); title('B298964');


sourcy = 'B162582';
source = ['/export/res/breast/' sourcy '/research_results_combined_kinetics/slice16/wash_out_vector.mat'];
load(source, 'wash_out_ROI_vec');
wo = wash_out_ROI_vec;
wash_out_vec_b = [wash_out_vec_b wash_out_ROI_vec];

sourcy = 'B162582';
source = ['/export/res/breast/' sourcy '/research_results_combined_kinetics/slice16/wash_in_vector_uptake.mat'];
load(source, 'wash_in_ROI_vec');
wi = wash_in_ROI_vec;
wash_in_vec_b = [wash_in_vec_b wi];

wo_mean = [];
wo_mean2 = [];
wi_vec = [];
for i=60:3:90 % walk through the bins
    ind = find(wi>i & wi<i+3);
    if numel(ind) ~= 0
    for j = 1:size(ind,1)
        tempo = ind(j);
        wo_mean2 = [wo_mean2 wo(tempo)];
    end  
    wi_vec(((i-60)/3)+1) = i+3;
    wo_mean(((i-60)/3)+1) = mean(wo_mean2);
    last = find(wob(((i-60)/3)+1,:)~=0);
    wob(((i-60)/3)+1,last:last+size(wo_mean2,1)) = wo_mean2; 
    wo_mean2 = [];
    end
end
% figure
% plot(wi_vec,wo_mean,'d'); xlabel('wi (degree)'), ylabel('wo (degree)'); title('B162582');

sourcy = 'B153241';
source = ['/export/res/breast/' sourcy '/research_results_combined_kinetics/slice5/wash_out_vector.mat'];
load(source, 'wash_out_ROI_vec');
wo = wash_out_ROI_vec;
wash_out_vec_b = [wash_out_vec_b wash_out_ROI_vec];

source = ['/export/res/breast/' sourcy '/research_results_combined_kinetics/slice5/wash_in_vector_uptake.mat'];
load(source, 'wash_in_ROI_vec');
wi = wash_in_ROI_vec;
wash_in_vec_b = [wash_in_vec_b wi];

wo_mean = [];
wo_mean2 = [];
wi_vec = [];
for i=60:3:90 % walk through the bins
    ind = find(wi>i & wi<i+3);
    if numel(ind) ~= 0
    for j = 1:size(ind,1)
        tempo = ind(j);
        wo_mean2 = [wo_mean2 wo(tempo)];
    end  
    wi_vec(((i-60)/3)+1) = i+3;
    wo_mean(((i-60)/3)+1) = mean(wo_mean2);
    last = find(wob(((i-60)/3)+1,:)~=0);
    wob(((i-60)/3)+1,last:last+size(wo_mean2,1)) = wo_mean2; 
    wo_mean2 = [];
    end
end
% figure
% plot(wi_vec,wo_mean,'d'); xlabel('wi (degree)'), ylabel('wo (degree)'); title('B153241');

sourcy = 'B49217';
source = ['/export/res/breast/' sourcy '/research_results_combined_kinetics/slice5/wash_out_vector.mat'];
load(source, 'wash_out_ROI_vec');
wo = wash_out_ROI_vec;
wash_out_vec_b = [wash_out_vec_b wash_out_ROI_vec];

source = ['/export/res/breast/' sourcy '/research_results_combined_kinetics/slice5/wash_in_vector_uptake.mat'];
load(source, 'wash_in_ROI_vec');
wi = wash_in_ROI_vec;
wash_in_vec_b = [wash_in_vec_b wi];

wo_mean = [];
wo_mean2 = [];
wi_vec = [];
for i=60:3:90 % walk through the bins
    ind = find(wi>i & wi<i+3);
    if numel(ind) ~= 0
    for j = 1:size(ind,1)
        tempo = ind(j);
        wo_mean2 = [wo_mean2 wo(tempo)];
    end  
    wi_vec(((i-60)/3)+1) = i+3;
    wo_mean(((i-60)/3)+1) = mean(wo_mean2);
        last = find(wob(((i-60)/3)+1,:)~=0);
    wob(((i-60)/3)+1,last:last+size(wo_mean2,1)) = wo_mean2; 
    wo_mean2 = [];
    end
end
% figure
% plot(wi_vec,wo_mean,'d'); xlabel('wi (degree)'), ylabel('wo (degree)'); title('B49217');
% 
% sourcy = 'B299986';
% source = ['/export/res/breast/' sourcy '/research_results_combined_kinetics/slice10/wash_out_vector.mat'];
% load(source, 'wash_out_ROI_vec');
% wo = wash_out_ROI_vec;
% wash_out_vec_b = [wash_out_vec_b wash_out_ROI_vec];
% 
% source = ['/export/res/breast/' sourcy '/research_results_combined_kinetics/slice10/wash_in_vector_uptake.mat'];
% load(source, 'wash_in_ROI_vec');
% wi = wash_in_ROI_vec;
% wash_in_vec_b = [wash_in_vec_b wi];
% 
% wo_mean = [];
% wo_mean2 = [];
% wi_vec = [];
% for i=60:3:90 % walk through the bins
%     ind = find(wi>i & wi<i+3);
%     if numel(ind) ~= 0
%     for j = 1:size(ind,1)
%         tempo = ind(j);
%         wo_mean2 = [wo_mean2 wo(tempo)];
%     end  
%     wi_vec(((i-60)/3)+1) = i+3;
%     wo_mean(((i-60)/3)+1) = mean(wo_mean2);
%         last = find(wob(((i-60)/3)+1,:)~=0);
%     wob(((i-60)/3)+1,last:last+size(wo_mean2,1)) = wo_mean2; 
%     wo_mean2 = [];
%     end
% end
% figure
% plot(wi_vec,wo_mean,'d'); xlabel('wi (degree)'), ylabel('wo (degree)'); title('B299986');

sourcy = 'N309953benign';
source = ['/export/res/breast/' sourcy '/research_results_combined_kinetics/slice9/wash_out_vector.mat'];
load(source, 'wash_out_ROI_vec');
wo = wash_out_ROI_vec;
wash_out_vec_b = [wash_out_vec_b wash_out_ROI_vec];

source = ['/export/res/breast/' sourcy '/research_results_combined_kinetics/slice9/wash_in_vector_uptake.mat'];
load(source, 'wash_in_ROI_vec');
wi = wash_in_ROI_vec;
wash_in_vec_b = [wash_in_vec_b wi];

wo_mean = [];
wo_mean2 = [];
wi_vec = [];
for i=60:3:90 % walk through the bins
    ind = find(wi>i & wi<i+3);
    if numel(ind) ~= 0
    for j = 1:size(ind,1)
        tempo = ind(j);
        wo_mean2 = [wo_mean2 wo(tempo)];
    end  
    wi_vec(((i-60)/3)+1) = i+3;
    wo_mean(((i-60)/3)+1) = mean(wo_mean2);
        last = find(wob(((i-60)/3)+1,:)~=0);
    wob(((i-60)/3)+1,last:last+size(wo_mean2,1)) = wo_mean2; 
    wo_mean2 = [];
    end
end
% figure
% plot(wi_vec,wo_mean,'d'); xlabel('wi (degree)'), ylabel('wo (degree)'); title('N309953');


size(wob)
for i=1:11
    wob(i,:)=mean(wob(i,:));
end
wi_vecb = 60:3:90;
% plot(wi_vecb, wob); title('dependence');


wi = wash_in_vec_b;
wo = wash_out_vec_b;
wo_mean = [];
wo_mean2 = [];
wi_vec = [];
index = 0;
for i=60:1:90 % walk through the bins
    index = index+1;
    ind = find(wi>i-1 & wi<i);
%     if numel(ind) ~= 0
    for j = 1:size(ind,2)
        tempo = ind(j);
        wo_mean2 = [wo_mean2 wo(tempo)];
    end  
%     wi_vec(index) = i+3;
    wo_mean(index) = mean(wo_mean2);
    size(wo_mean2)
%         last = find(wob(((i-60)/3)+1,:)~=0);
%     wob(index,last:last+size(wo_mean2,1)) = wo_mean2; 
    wo_mean2 = [];
%     end
end
% figure
% plot(wi_vec,wo_mean,'d'); xlabel('wi (degree)'), ylabel('wo (degree)'); title('N309953');


% size(wob)
% for i=1:11
%     wob(i,:)=mean(wob(i,:));
% end
% wi_vecb = 60:1:90;
% figure; plot(wi_vecb, wo_mean,'d'); title('dependence');


% sourcy = 'B162582';
% source = ['/export/res/breast/' sourcy '/research_results_combined_kinetics/slice16/wash_out_vector.mat'];
% load(source, 'wash_out_ROI_vec');
% wash_out_vec_b = [wash_out_vec_b wash_out_ROI_vec];
% 
% sourcy = 'B153241';
% source = ['/export/res/breast/' sourcy '/research_results_combined_kinetics/slice5/wash_out_vector.mat'];
% load(source, 'wash_out_ROI_vec');
% wash_out_vec_b = [wash_out_vec_b wash_out_ROI_vec];
% 
% sourcy = 'B49217';
% source = ['/export/res/breast/' sourcy '/research_results_combined_kinetics/slice5/wash_out_vector.mat'];
% load(source, 'wash_out_ROI_vec');
% wash_out_vec_b = [wash_out_vec_b wash_out_ROI_vec];
% 
% sourcy = 'B299986';
% source = ['/export/res/breast/' sourcy '/research_results_combined_kinetics/slice10/wash_out_vector.mat'];
% load(source, 'wash_out_ROI_vec');
% wash_out_vec_b = [wash_out_vec_b wash_out_ROI_vec];
% 
% sourcy = 'N309953benign';
% source = ['/export/res/breast/' sourcy '/research_results_combined_kinetics/slice9/wash_out_vector.mat'];
% load(source, 'wash_out_ROI_vec');
% wash_out_vec_b = [wash_out_vec_b wash_out_ROI_vec];
% 
% 
sourcy = 'M298821';
source = ['/export/res/breast/' sourcy '/research_results_combined_kinetics/slice29/wash_out_vector.mat'];
load(source, 'wash_out_ROI_vec');
wo = wash_out_ROI_vec;
wash_out_vec_m = [wash_out_vec_m wash_out_ROI_vec];

source = ['/export/res/breast/' sourcy '/research_results_combined_kinetics/slice29/wash_in_vector_uptake.mat'];
load(source, 'wash_in_ROI_vec');
wi = wash_in_ROI_vec;
wash_in_vec_m = [wash_in_vec_m wi];

wo_mean = [];
wo_mean2 = [];
wi_vec = [];
for i=60:3:90 % walk through the bins
    ind = find(wi>i & wi<i+3);
    if numel(ind) ~= 0
    for j = 1:size(ind,1)
        tempo = ind(j);
        wo_mean2 = [wo_mean2 wo(tempo)];
    end  
    wi_vec(((i-60)/3)+1) = i+3;
    wo_mean(((i-60)/3)+1) = mean(wo_mean2);
    wo_mean2 = [];
    end
end
% figure
% plot(wi_vec,wo_mean,'d'); xlabel('wi (degree)'), ylabel('wo (degree)'); title('M298821');


sourcy = 'M301431';
source = ['/export/res/breast/' sourcy '/research_results_combined_kinetics/slice19/wash_out_vector.mat'];
load(source, 'wash_out_ROI_vec');
wo = wash_out_ROI_vec;
wash_out_vec_m = [wash_out_vec_m wash_out_ROI_vec];

source = ['/export/res/breast/' sourcy '/research_results_combined_kinetics/slice19/wash_in_vector_uptake.mat'];
load(source, 'wash_in_ROI_vec');
wi = wash_in_ROI_vec;
wash_in_vec_m = [wash_in_vec_m wi];

sourcy = 'M179579';
source = ['/export/res/breast/' sourcy '/research_results_combined_kinetics/slice11/wash_out_vector.mat'];
load(source, 'wash_out_ROI_vec');
wo = wash_out_ROI_vec;
wash_out_vec_m = [wash_out_vec_m wash_out_ROI_vec];

source = ['/export/res/breast/' sourcy '/research_results_combined_kinetics/slice11/wash_in_vector_uptake.mat'];
load(source, 'wash_in_ROI_vec');
wi = wash_in_ROI_vec;
wash_in_vec_m = [wash_in_vec_m wi];

sourcy = 'M7803';
source = ['/export/res/breast/' sourcy '/research_results_combined_kinetics/slice21/wash_out_vector.mat'];
load(source, 'wash_out_ROI_vec');
wo = wash_out_ROI_vec;
wash_out_vec_m = [wash_out_vec_m wash_out_ROI_vec];

source = ['/export/res/breast/' sourcy '/research_results_combined_kinetics/slice21/wash_in_vector_uptake.mat'];
load(source, 'wash_in_ROI_vec');
wi = wash_in_ROI_vec;
wash_in_vec_m = [wash_in_vec_m wi];

sourcy = 'M273920';
source = ['/export/res/breast/' sourcy '/research_results_combined_kinetics/slice8/wash_out_vector.mat'];
load(source, 'wash_out_ROI_vec');
wo = wash_out_ROI_vec;
wash_out_vec_m = [wash_out_vec_m wash_out_ROI_vec];

source = ['/export/res/breast/' sourcy '/research_results_combined_kinetics/slice8/wash_in_vector_uptake.mat'];
load(source, 'wash_in_ROI_vec');
wi = wash_in_ROI_vec;
wash_in_vec_m = [wash_in_vec_m wi];

sourcy = 'N84432';
source = ['/export/res/breast/' sourcy '/research_results_combined_kinetics/slice8/wash_out_vector.mat'];
load(source, 'wash_out_ROI_vec');
wo = wash_out_ROI_vec;
wash_out_vec_m = [wash_out_vec_m wash_out_ROI_vec];

source = ['/export/res/breast/' sourcy '/research_results_combined_kinetics/slice8/wash_in_vector_uptake.mat'];
load(source, 'wash_in_ROI_vec');
wi = wash_in_ROI_vec;
wash_in_vec_m = [wash_in_vec_m wi];

sourcy = 'M297672';
source = ['/export/res/breast/' sourcy '/research_results_combined_kinetics/slice12/wash_out_vector.mat'];
load(source, 'wash_out_ROI_vec');
wo = wash_out_ROI_vec;
wash_out_vec_m = [wash_out_vec_m wash_out_ROI_vec];

source = ['/export/res/breast/' sourcy '/research_results_combined_kinetics/slice12/wash_in_vector_uptake.mat'];
load(source, 'wash_in_ROI_vec');
wi = wash_in_ROI_vec;
wash_in_vec_m = [wash_in_vec_m wi];

sourcy = 'B49217b';
source = ['/export/res/breast/' sourcy '/research_results_combined_kinetics/slice7/wash_out_vector.mat'];
load(source, 'wash_out_ROI_vec');
wo = wash_out_ROI_vec;
wash_out_vec_m = [wash_out_vec_m wash_out_ROI_vec];

source = ['/export/res/breast/' sourcy '/research_results_combined_kinetics/slice7/wash_in_vector_uptake.mat'];
load(source, 'wash_in_ROI_vec');
wi = wash_in_ROI_vec;
wash_in_vec_m = [wash_in_vec_m wi];

sourcy = 'M303424';
source = ['/export/res/breast/' sourcy '/research_results_combined_kinetics/slice16/wash_out_vector.mat'];
load(source, 'wash_out_ROI_vec');
wo = wash_out_ROI_vec;
wash_out_vec_m = [wash_out_vec_m wash_out_ROI_vec];

source = ['/export/res/breast/' sourcy '/research_results_combined_kinetics/slice16/wash_in_vector_uptake.mat'];
load(source, 'wash_in_ROI_vec');
wi = wash_in_ROI_vec;
wash_in_vec_m = [wash_in_vec_m wi];

sourcy = 'N309953';
source = ['/export/res/breast/' sourcy '/research_results_combined_kinetics/slice11/wash_out_vector.mat'];
load(source, 'wash_out_ROI_vec');
wo = wash_out_ROI_vec;
wash_out_vec_m = [wash_out_vec_m wash_out_ROI_vec];

source = ['/export/res/breast/' sourcy '/research_results_combined_kinetics/slice11/wash_in_vector_uptake.mat'];
load(source, 'wash_in_ROI_vec');
wi = wash_in_ROI_vec;
wash_in_vec_m = [wash_in_vec_m wi];

% 
% sourcy = 'M301431';
% source = ['/export/res/breast/' sourcy '/research_results_combined_kinetics/slice19/wash_out_vector.mat'];
% load(source, 'wash_out_ROI_vec');
% wash_out_vec_m = [wash_out_vec_m wash_out_ROI_vec];
% 
% sourcy = 'M179579';
% source = ['/export/res/breast/' sourcy '/research_results_combined_kinetics/slice11/wash_out_vector.mat'];
% load(source, 'wash_out_ROI_vec');
% wash_out_vec_m = [wash_out_vec_m wash_out_ROI_vec];
% 
% sourcy = 'M7803';
% source = ['/export/res/breast/' sourcy '/research_results_combined_kinetics/slice21/wash_out_vector.mat'];
% load(source, 'wash_out_ROI_vec');
% wash_out_vec_m = [wash_out_vec_m wash_out_ROI_vec];
% 
% sourcy = 'M273920';
% source = ['/export/res/breast/' sourcy '/research_results_combined_kinetics/slice8/wash_out_vector.mat'];
% load(source, 'wash_out_ROI_vec');
% wash_out_vec_m = [wash_out_vec_m wash_out_ROI_vec];
% 
% sourcy = 'N84432';
% source = ['/export/res/breast/' sourcy '/research_results_combined_kinetics/slice8/wash_out_vector.mat'];
% load(source, 'wash_out_ROI_vec');
% wash_out_vec_m = [wash_out_vec_m wash_out_ROI_vec];
% 
% sourcy = 'M297672';
% source = ['/export/res/breast/' sourcy '/research_results_combined_kinetics/slice12/wash_out_vector.mat'];
% load(source, 'wash_out_ROI_vec');
% wash_out_vec_m = [wash_out_vec_m wash_out_ROI_vec];
% 
% sourcy = 'B49217b';
% source = ['/export/res/breast/' sourcy '/research_results_combined_kinetics/slice7/wash_out_vector.mat'];
% load(source, 'wash_out_ROI_vec');
% wash_out_vec_m = [wash_out_vec_m wash_out_ROI_vec];
% 
% sourcy = 'M303424';
% source = ['/export/res/breast/' sourcy '/research_results_combined_kinetics/slice16/wash_out_vector.mat'];
% load(source, 'wash_out_ROI_vec');
% wash_out_vec_m = [wash_out_vec_m wash_out_ROI_vec];
% 
% sourcy = 'N309953';
% source = ['/export/res/breast/' sourcy '/research_results_combined_kinetics/slice11/wash_out_vector.mat'];
% load(source, 'wash_out_ROI_vec');
% wash_out_vec_m = [wash_out_vec_m wash_out_ROI_vec];
% 
% temp = wash_out_vec_m(wash_out_vec_m~=0);
% mal = temp;
% mean_malignant = mean(temp)
% std_malignant = std(temp)
% figure; hist(temp,50); title('malignant');
% % figure; histfit(temp,50); title('histogram incl. fit');
% 
% wash_out_mean = mean(wash_out_vec_m(wash_out_vec_m~=0))
% wash_out_std = std(wash_out_vec_m(wash_out_vec_m~=0))
% 
% x = -90:3:90;
% 
% d = (x-wash_out_mean)/wash_out_std;
% % gauss = 18474*(1/sqrt(2*pi))*exp(-d.^2/2);
% gauss = (1/sqrt(2*pi))*exp(-d.^2/2)/wash_out_std;
% figure; plot(x,gauss);
% hold on
% 
% h_m = hist(temp,x);
% hn_m = h_m/(length(temp)*3); % 3 is the bin width
% bar(x,hn_m);
% hold off;
% 
% figure
% temp = wash_out_vec_b(wash_out_vec_b~=0);
% ben = temp;
% mean_benign = mean(temp)
% std_benign = std(temp)
% h_b = hist(temp,x);
% hn_b = h_b/(length(temp)*3); 
% bar(x,hn_b);
% 
% 
% % figure; hist(wash_out_vec_b(wash_out_vec_b~=0),50); title('benign');
% 
% % figure; hist(temp,x); title('malignant distribution');
% % plot(x,gauss);
% % hold off
% 
% 
% 
% % larger=numel(wash_out_vec(wash_out_vec>0))
% % smaller=numel(wash_out_vec(wash_out_vec<0))
% 
% % wash_out_vec = [wash_out_vec_b wash_out_vec_m];
% % figure; hist(wash_out_vec(wash_out_vec~=0),50); title('all lesions');
