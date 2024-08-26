function eyeana(file,seq,run_ind)
%eyeana(file,seq,run_ind)
if ~exist('run_ind','var')
    run_ind = 2:11;
end

load(file);

hfactor = (atan(23/78)/pi*180)/(215-42);
vfactor = (atan(18/78)/pi*180)/(193-46);
centerX = 125;
centerY = 122;

fix_percent = zeros(1,length(record_arr));
var_fixx = zeros(1,length(record_arr));
var_fixy = zeros(1,length(record_arr));
var_x = zeros(1,length(record_arr));
var_y = zeros(1,length(record_arr));

for i=1:length(record_arr)  
  data = record_arr{i};
  x_pos =  (data(11,:)-centerX)*hfactor;
  y_pos =  (data(12,:)-centerY)*vfactor;
 % x_pos = detrend(x_pos);
 % y_pos = detrend(y_pos);
  tmp = find(abs(x_pos)<5 & abs(y_pos)<5);
  var_fixx(i) = std(x_pos(tmp));
  var_fixy(i) = std(y_pos(tmp));
  var_x(i) = std(x_pos);
  var_y(i) = std(y_pos);
  
  fix_percent(i) = length(tmp)/length(x_pos);
 % figure;
 % subplot(2,1,1);plot(x_pos,'k-');
 % subplot(2,1,2);
 % plot(y_pos,'k-');
end
fprintf('%3.2f ',fix_percent);fprintf('\n');
fprintf('%3.2f ',var_fixx);fprintf('\n');
fprintf('%3.2f ',var_fixy);fprintf('\n');
fprintf('%3.2f ',var_x);fprintf('\n');
fprintf('%3.2f ',var_y);fprintf('\n');

if ~exist('seq','var')
    return;
end
seq = load(seq);
std_x = zeros(length(seq(:))/4,4);
std_y = zeros(length(seq(:))/4,4);
percent_fix = zeros(length(seq(:))/4,4);

icoh = zeros(1,4);

x_posc = cell(1,4);
y_posc = cell(1,4);

for i=1:length(run_ind)
    
  data = record_arr{run_ind(i)};

  xd = data(9,:);
  ind = find(xd>4);
  dind = diff(ind);
  nTR = length(find(dind>1))+1;
  disp(nTR);
  ind_TR(1) = ind(1);
  ind_TR(2:nTR) = ind(find(dind>1)+1);

  
  for itrl=1:8
      coh = seq(i,itrl);
      i1 = ind_TR(10+(itrl-1)*15);
      if itrl <8
          i2 = ind_TR(10+itrl*15)-1;
      else
          i2 = i1+119;
      end
      x_pos = (data(11,:)-centerX)*hfactor;
      y_pos = (data(12,:)-centerY)*vfactor;
      mask = abs(x_pos)<2 & abs(y_pos)<2;
      x_mn = mean(x_pos(mask));
      y_mn = mean(y_pos(mask));
      
      x_tmp = (data(11,i1:i2)-centerX)*hfactor;
      y_tmp = (data(12,i1:i2)-centerY)*vfactor;
      mask = abs(x_tmp)>2 | abs(y_tmp)>2;
      
      icoh(coh) = icoh(coh)+1;
      percent_fix(icoh(coh),coh) = 1-length(find(mask))/length(x_tmp);
      x_tmp(mask) = [];
      y_tmp(mask) = [];
      
      std_x(icoh(coh),coh) = sqrt(mean((x_tmp-x_mn).^2));
      std_y(icoh(coh),coh) = sqrt(mean((y_tmp-y_mn).^2));
      x_posc{coh}(end+1:end+length(x_tmp)) =  x_tmp-x_mn;
      y_posc{coh}(end+1:end+length(y_tmp)) =  y_tmp-y_mn;
      
  end
      
end
disp('fixation deviation in x');
disp(mean(std_x,1));
disp(std(std_x,0,1));

disp('fixation deviation in y');
disp(mean(std_y,1));
disp(std(std_y,0,1));
disp('F test for difference between the four conditions');
disp(ftest_matrx(std_x));
disp(ftest_matrx(std_y))

disp('percent fixation');
disp(mean(percent_fix,1));
disp(std(percent_fix,0,1));
disp(ftest_matrx(percent_fix));

for i=1:4
  figure;
  subplot(2,1,1);plot(x_posc{i},'k-');
  subplot(2,1,2);
  plot(y_posc{i},'k-');
  disp([std(x_posc{i}),std(y_posc{i})]);
end
