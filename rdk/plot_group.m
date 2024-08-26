function plot_group
global pk;
pk = 4:8;

root = '/home/zong/rdk';
flist = {'Zong/rdk_12_01_08/afni','Peng/RDK_12_16_2008/afni','Hua/rdk_12_04_08/afni','KELIRIS004_GeorgeSister/rdk/afni','KELIRIS003/afni','KELIRIS005/rdk/afni'};
%dlist = {'results_xcoh','results_xcoh','results_xcoh','thr1.92','results_xcoh_thr1.92','results_xcoh'};
%dlist = {'.','.','.','.'};
for i=1:6
  dlist{i} = fullfile(root,flist{i});
end

%plot_V5_anpos(dlist(1:3),[1,1,1]);
%plot_V5(dlist,[1,1,1,0]);
%plot_h_V1to4(dlist,[1,1,1,0]);
%plot_ts_V1to4(dlist,[1,1,1,0]);
stat_test(dlist(1:4),[1,1,1,0],[0,0,0,0]);
%stat_test(dlist(1:3),[1,1,1],[1,1,1]);
%stat_test_glm(dlist(1:4),[1,1,1,0],[0,0,0,0]);
%plot_V1to5(dlist,[1,1,1,0]);
%plot_V1to4ecc(dlist,[1,1,1,0]);
%plot_glm_rois(dlist(1:4),[1,1,1,0]);
function [p,stat]=do_anova1(h_arr,coh4)
 d=[];
 group = [];
 ns = size(h_arr,2);
  for i=1:size(h_arr,1)
      if i==1
          d = [d,h_arr(i,coh4>0)];
          group = [group,ones(1,length(find(coh4>0)))*i];
      else
          d = [d,h_arr(i,:)];
          group = [group,ones(1,ns)*i];
      end
  end
  
  [p,table,stat]=anova1(d,group,'off');
  
function p=do_anova2(h_arr,coh4)
 d=[];
 sid = [];
 coh = [];
 ns = size(h_arr,2);
  for i=1:size(h_arr,1)
      for j=1:ns
        if i~=1 || coh4(j) >0
          if iscell(h_arr)
            n = length(h_arr{i,j});  
            d = [d,h_arr{i,j}];
            sid = [sid,ones(1,n)*j];
            coh = [coh,ones(1,n)*i];
          else
            d = [d,h_arr(i,j)];
            sid = [sid,j];
            coh = [coh,i];
          end
        end
      end
  end
  
  [p,table,stat]=anovan(d,{sid,coh},'model',2,'sstype',3,'varname',{'sid';'coh'},'display','off');
  
  figure;c=multcompare(stat,'dimension',2);
  %disp(c);
  

function plot_glm_rois(dlist,coh4)


[a,astd,h,eh,h_arr]=groupAverage(dlist(1:4),'ts_V5glm.mat',false,coh4);
%ST 1,2,4
figure;

set(gcf,'Units','inch','Position',[2,3,6.2,8]);
subplot(4,2,1);plot_hght(h,eh);
ylim([0,1.3]);text(2,0.5,'hMT+','FontSize',16);
subplot(4,2,2);
[a,astd,h,eh,h_arr]=groupAverage(dlist([1,3,4]),'ts_STglm.mat',false,coh4([1,3,4]));
plot_hght(h,eh);
ylim([0,1.1]);text(2,0.5,'ST','FontSize',16);
subplot(4,2,3);
[a,astd,h,eh,h_arr]=groupAverage(dlist,'ts_V3Aglm.mat',false,coh4);
plot_hght(h,eh);
ylim([0,1.2]);text(2,0.5,'Cun','FontSize',16);
subplot(4,2,4);
[a,astd,h,eh,h_arr]=groupAverage(dlist([1,2,4]),'ts_LSglm.mat',false,coh4([1,2,4]));
plot_hght(h,eh);
ylim([0,1]);text(2,0.5,'pLS','FontSize',16);

[a,astd,h,eh,h_arr]=groupAverage(dlist([1,2]),'ts_IPSglm.mat',false,coh4([1,2]));
subplot(4,2,5); plot_hght(h,eh);
ylim([0,1.6]);text(2,0.5,'IPS','FontSize',16);
[a,astd,h1]=get_data(dlist{1},'ts_CCglm.mat');
subplot(4,2,6); plot_hght(h1',zeros(size(h1')));
ylim([0,1.3]);text(2,0.5,'Cing','FontSize',16);
%xlabel('Coherence (%)');
text(2,-0.5,'Coherence (%)','FontSize',16);
[a,astd,h2]=get_data(dlist{2},'ts_PreCglm.mat');
subplot(4,2,7); plot_hght(h2',zeros(size(h2')));
ylim([0,1.1]);text(2,0.5,'PreC','FontSize',16);

set(gca,'FontSize',16);
text(2,-0.5,'Coherence (%)','FontSize',16);

%xlabel('Coherence (%)');
ylabel('Height (%)');

function plot_V5_anpos(dlist,coh4)

figure;

  [a,astd,h,eh]=groupAverage(dlist(1:3),'ts_V5u_a.mat',false,coh4);
  [a2,astd2,h2,eh2]=groupAverage(dlist(1:3),'ts_V5u_p.mat',false,coh4);
  
plot_hght_bar(h,eh,h2,eh2);
   p=calc_p(h,eh,h2,eh2,[3,3,3,3]);
   disp(p);
   
   ylim([0,2]);
    xlabel('Coherence (%)');

    ylabel('Height (%)');

set(gca,'FontSize',16);
text(2,0.2,'hMT+','FontSize',16);
set(gca,'XTickLabel',{'12.5','25','50','100'});

for i=1:4
    figure;
    errorbar(a(:,i),astd(:,i));
    hold on;
    errorbar(a2(:,i),astd(:,i),'r');
    title(num2str(i));
end

function plot_h_V1to4(dlist,coh4)

[a,astd,h1,eh1,h_arr1]=groupAverage(dlist(1:4),'ts_V1u.mat',false,[1,1,1,0]);
[a,astd,h2,eh2,h_arr2]=groupAverage(dlist(1:4),'ts_V2u.mat',false,[1,1,1,0]);
[a,astd,h3,eh3,h_arr3]=groupAverage(dlist(1:4),'ts_V3u.mat',false,[1,1,1,0]);
[a,astd,h4,eh4,h_arr4]=groupAverage(dlist(1:4),'ts_V4u.mat',false,[1,1,1,0]);

figure;

set(gcf,'Units','inch','Position',[3,4,5.8,4.2]);
subplot(2,2,1);
plot_hght(h1,eh1,h_arr1,coh4);
ylim([-1.2,0.4]);
title('V1');
set(gca,'XTickLabel','');
set(gca,'YTick',-1.2:0.4:1);

subplot(2,2,2);
plot_hght(h2,eh2,h_arr2,coh4);
title('V2');
set(gca,'YTick',-1.2:0.4:1);
set(gca,'XTickLabel','');

subplot(2,2,3);
plot_hght(h3,eh3,h_arr3,coh4);
ylabel('Height (%)');
xlabel('Coherence (%)');
title('V3');
set(gca,'YTick',-0.9:0.3:1);

subplot(2,2,4);
plot_hght(h3,eh4,h_arr4,coh4);
xlabel('Coherence (%)');
title('V4');
set(gca,'YTick',-0.9:0.3:1)

function stat_test(dlist,coh4,coh4b)

for i=1:5
   [a,astd,h,eh,h_arr,h_arr2]=groupAverage(dlist,sprintf('ts_V%du.mat',i),false,coh4);
   % h_arr: [ncoh,nsub];
   % h_arr2: cell{ncoh,nsub};
   show_p(h_arr,h_arr2,coh4b,sprintf('V%d',i));
   
end

function stat_test_glm(dlist,coh4,coh4b)

% coh4b: whether to include 12.5% in the statistical test. 

[a,astd,h,eh,h_arr,h_arr2]=groupAverage(dlist(1:4),'ts_V5glm.mat',false,coh4);
show_p(h_arr,h_arr2,coh4b,'V5glm');

[a,astd,h,eh,h_arr,h_arr2]=groupAverage(dlist([1,3,4]),'ts_STglm.mat',false,coh4([1,3,4]));
show_p(h_arr,h_arr2,coh4b([1,3,4]),'STglm');

[a,astd,h,eh,h_arr,h_arr2]=groupAverage(dlist,'ts_V3Aglm.mat',false,coh4);
show_p(h_arr,h_arr2,coh4b,'V3Aglm');

[a,astd,h,eh,h_arr,h_arr2]=groupAverage(dlist([1,2,4]),'ts_LSglm.mat',false,coh4([1,2,4]));
show_p(h_arr,h_arr2,coh4b([1,2,4]),'LSglm');

[a,astd,h,eh,h_arr,h_arr2]=groupAverage(dlist([1,2]),'ts_IPSglm.mat',false,coh4([1,2]));
show_p(h_arr,h_arr2,coh4b([1,2]),'IPSglm');

[a,astd,h1,h_arr2]=get_data(dlist{1},'ts_CCglm.mat');
[p,stat]=do_anova1_ss(h_arr2,coh4b(1));
fprintf('CCglm: P value: %5.4f\n\n',p);
[a,astd,h2,h_arr2]=get_data(dlist{2},'ts_PreCglm.mat');
[p,stat]=do_anova1_ss(h_arr2,coh4b(2));
fprintf('PreCglm: P value: %5.4f\n\n',p);

   
function show_p(h_arr,h_arr2,coh4,roi)
   
   [p,stat]=do_anova1(h_arr,coh4);
   fprintf('%s: P values: %5.4f\n',roi,p);
   figure;multcompare(stat);
   title(roi);
   
   fprintf('Single subject: ');
   for j=1:length(coh4)
     [p,stat] = do_anova1_ss(h_arr2(:,j),coh4(j));
     fprintf('%5.4f ;',p);
   end 
   fprintf('\n\n');
   %sig_diff_zero(h,eh,coh4); 
        

function [p,stat]=do_anova1_ss(h_arr,coh4)

d = [];
group = [];
    for i=1:length(h_arr)
        if coh4 || i>1 
          d = [d,h_arr{i}];
          group = [group,i*ones(1,length(h_arr{i}))];
        end  
    end
   [p,table,stat]= anova1(d,group,'off');
    
    
        
function plot_V1to5(dlist,coh4)

for i=1:5
   [a{i},astd{i},h{i},eh{i},h_arr{i}]=groupAverage(dlist(1:4),sprintf('ts_V%du.mat',i),false,coh4);
end
%calc_p(h(1),eh(1),h(3),eh(3),4);

figure;
set(gcf,'Units','inch','Position',[2,3,6.2,10]);
yl = [-0.7,0.5;-0.5,0.5;0,0.4;0,0.4;0,1];
yl_ts=[-1,0.5;-0.5,0.5;-0.2,0.4;-0.2,0.6;-0.3,1];
lb = {'V1','V2','V3','V4','hMT+'};
for i=1:5
subplot(5,2,(i-1)*2+1);
plot_ts(a{i},astd{i},yl_ts(i,:),lb{i});
subplot(5,2,(i-1)*2+2);
plot_hght(h{i},eh{i},h_arr{i},coh4);
text(2,0.5,lb{i},'FontSize',16);
ylim(yl(i,:));
end
subplot(5,2,10);ylabel('Height (%)');ylim([0,1.1]);
subplot(5,2,10);
text(2,0.5,'Coherence (%)','FontSize',16);
subplot(5,2,9);ylabel('Signal Change (%)');
subplot(5,2,9);%
text(10,0.5,'Time (2 s)','FontSize',16);

subplot(5,2,6);set(gca,'YTick',-0.2:0.2:0.6);
subplot(5,2,8);set(gca,'YTick',-0.2:0.2:0.6);


subplot(5,2,1);
a=findobj(gca,'Type','Line');
legend(gca,a([2,1]),{'50%','100%'},'FontSize',16);

ah=axes('Position',get(gca,'Position'),'Visible','off');
legend(ah,a([4,3]),{'12.5%','25%'},'FontSize',16)

%legend boxoff;

function plot_V1to4ecc(dlist,coh4)

figure;

for i=1:4
  [a{i},astd{i},h_sm,eh_sm]=groupAverage(dlist(1:4),sprintf('ts_V%du_se.mat',i),false,coh4);
  [a2{i},astd2{i},h_lg,eh_lg]=groupAverage(dlist(1:4),sprintf('ts_V%du_le.mat',i),false,coh4);
  subplot(2,2,i);
plot_hght_bar(h_sm,eh_sm,h_lg,eh_lg);
p=calc_p(h_sm,eh_sm,h_lg,eh_lg,[3,4,4,4]);
disp(i);
disp(p);
if i==3 || i==4
    ylim([0,0.4]);
    xlabel('Coherence (%)');
end
if i==3
    ylabel('Signal Change (%)');
end
set(gca,'FontSize',16);
text(2,0.2,sprintf('V%d',i),'FontSize',16);
set(gca,'XTickLabel',{'12.5','25','50','100'});

end

for i=1:4
    figure;
    errorbar(a{i}(:,2),astd{i}(:,2));
    hold on;
    errorbar(a2{i}(:,2),astd{i}(:,2),'r');
    title(sprintf('V%d',i));
end


function plot_ts_V1to4(dlist,coh4)

[a1,astd1]=groupAverage(dlist(1:4),'ts_V1u.mat',false,coh4);
[a2,astd2]=groupAverage(dlist(1:4),'ts_V2u.mat',false,coh4);
[a3,astd3]=groupAverage(dlist(1:4),'ts_V3u.mat',false,coh4);
[a4,astd4]=groupAverage(dlist(1:4),'ts_V4u.mat',false,coh4);

figure;

set(gcf,'Units','inch','Position',[3,4,5.8,4.2]);
subplot(2,2,1);
plot_ts(a1,astd1,[-1,0.5],'V1');
set(gca,'YTick',-1:0.5:2);
set(gca,'XTickLabel','');
%title('V1');

subplot(2,2,2);
plot_ts(a2,astd2,[-1,0.5],'V2');
set(gca,'YTick',-1:0.5:2);
%title('V2');
set(gca,'XTickLabel','');
subplot(2,2,3);
plot_ts(a3,astd3,[-0.4,0.6],'V3');
set(gca,'YTick',-1:0.2:2);
ylabel('Signal Change (%)');
xlabel('Time (2 s)');
%title('V3');
subplot(2,2,4);
plot_ts(a4,astd4,[-0.4,0.6],'V4');
set(gca,'YTick',-1:0.2:2);
xlabel('Time (2 s)');
%title('V4');
legend({'12.5%','25%','50%','100%'},'FontSize',14);
%xlabel('t (2s)');
%ylabel('Signal Change (%)');
%legend('Group','Individual');

function plot_hght(h,eh,h_arr,coh4)

%figure;

nh = size(h,2);
for i=1:nh
   % errorbar([0.125,0.25,0.5,1]*100,h(:,i),eh(:,i),'ko-','MarkerSize',6,'LineWidth',2);
     hd = barweb(h,eh,0.7);
   set(hd.bars(1),'FaceColor',[0.5,0.5,0.5]); 
   hold on;
end

if exist('h_arr','var') && ~isempty(h_arr)
    for i=1:size(h_arr,2)
         if coh4(i)
   %        plot([0.125,0.25,0.5,1]*100,h_arr(:,i),'k>-','MarkerSize',7);
         else
    %       plot([0.25,0.5,1]*100,h_arr(2:4,i),'k>-','MarkerSize',7);
         end
     end
end
set(gca,'FontSize',16);
xlim([0.5,4.4]);
yrange =min_max([h+eh,h-eh]);
ylim(yrange+[-0.1,0.1]);

set(gca,'TickLength',[0.03,0.03])
set(gca,'XTick',1:1:5);
set(gca,'XTickLabel',{'12.5','25','50','100'});
set(gca,'YTick',-1:0.5:2);
set(gca,'YMinorTick','off');

function plot_hght_bar(h,eh,h2,eh2)

%figure;

    data = [h,h2];
    e = [eh,eh2];
   hd = barweb(data,e);
   set(hd.bars(1),'FaceColor',[0.5,0.5,0.5]);
   set(hd.bars(2),'FaceColor',[0.2,0.2,0.2]);
  
set(gca,'FontSize',16);

yrange =min_max([h+eh,h-eh]);
ylim(yrange+[-0.1,0.1]);

set(gca,'TickLength',[0.03,0.03])


function plot_ts(a,a_std,yrange,mask)

symb = {'k->','k-','k-o','k-s'};
dur = size(a,1);
for i=1:4
  %  errorbar(0:dur-1,a(:,i),a_std(:,i),symb{i},'LineWidth',1);
    plot(0:dur-1,a(:,i),symb{i},'LineWidth',1,'MarkerSize',3.5);
    hold on;
end
set(gca,'FontSize',16);
%yrange =min_max(cat(1,a+a_std,a-a_std));

xlim([0,dur]);
ylim(yrange);
rectangle('Position',[0,yrange(1),5,0.1*diff(yrange)],'FaceColor',[0.5,0.5,0.5]);
set(gca,'TickLength',[0.03,0.03]);
set(gca,'XTick',0:5:25);
set(gca,'Box','off');
text(15,0.2,mask,'FontSize',16);

function p=calc_p(h1,eh1,h2,eh2,n)
p = zeros(1,length(h1));
if numel(n) ==1
    n = n*ones(1,length(h1));
end
for i=1:length(h1)
    p(i) = ftest_mnstd([h1(i),h2(i)],[eh1(i),eh2(i)],n(i));
end
    
function sig_diff_zero(h,e,coh4)
% calculate whether h is signifcantly different from zero
p=zeros(size(h));
dof = length(coh4)-1;
for i=1:length(h)
    t = h(i)./e(i);
    if i==1
     dof = length(find(coh4>0));
    end
    p(i) = tTest(dof,t);
end

disp(p');
    

function [a,a_std,h,h_arr]=get_data(d,fname)
global pk;
var = 'ts_mn';
evar = 'ts_sem';
       
     if exist(fullfile(d,fname),'file') 
        load(fullfile(d,fname),var,evar,'n','ts');
     else
        error('File not found');
     end
h_arr = cell(4,1);

for i=1:4
    if length(ts) == 4
     b = 100*ts_detrend(ts{i},[1,16:20],0);
     h_arr{i} = mean(b(pk,:),1);
    elseif i>1
     b = 100*ts_detrend(ts{i-1},[1,16:20],0);
     h_arr{i} = mean(b(pk,:),1);
    else
    end
       
       
end
    
     a = eval(var);
 a_std = eval(evar);
 %arr = ts_detrend(arr,[1:4,21:24],1);
 a = 100*ts_detrend(a,[1,16:20],0);
 h = squeeze(mean(a(pk,:),1));
 
 
     

function [a,a_std,h,eh,h_arr,h_arr2] = groupAverage(dlist,fname,scale,coh4)
 global pk;
 coh4 = logical(coh4);
 nsub = length(dlist);
 var = 'ts_mn';
 
 for i=1:nsub
       
     if exist(fullfile(dlist{i},fname),'file') 
        load(fullfile(dlist{i},fname),var,'ts','n');
     else
        error(['File ', fullfile(dlist{i},fname),' not found']);
     end
    if i==1
        arr = zeros([size(eval(var)),nsub]);
        arr2 = cell(4,nsub);
        h_arr2 = cell(4,nsub);
    end
   
       if ~coh4(i)
         arr(:,2:4,i) = eval(var);
         for j=1:3
           arr2{j+1,i} = 100*ts_detrend(ts{j},[1,16:20],0);
           h_arr2{j+1,i} = mean(arr2{j+1,i}(pk,:),1);
         end
       else
          arr(:,:,i) = eval(var);
          for j=1:4
           arr2{j,i} = 100*ts_detrend(ts{j},[1,16:20],0);
           h_arr2{j,i} = mean(arr2{j,i}(pk,:),1);
         end
       end
 end

 arr = 100*ts_detrend(arr,[1,16:20],0);

 if scale
     h_arr = mean(arr(pk,:,:),1)-mean(arr(11:12,:,:),1);
     h_arr = abs(h_arr);
     arr = arr./repmat(h_arr(1,end,:),[size(arr,1),size(arr,2),1]); 
 end
     [h,eh,h_arr] = calc_hght(arr,coh4);
     
 a = mean(arr,3);
 a_std = std(arr,0,3)/sqrt(nsub);
 a(:,1) = mean(arr(:,1,coh4),3);
 a_std(:,1) = std(arr(:,1,coh4),0,3)/sqrt(length(find(coh4)));
 
 
function [h,eh,h_arr]=calc_hght(ts,coh4)
global pk;
sz = size(ts);
   nsub = sz(3);
  % h_arr = squeeze(mean(ts(pk,:,:),1)-mean(ts(11:12,:,:),1));
   
    h_arr = squeeze(mean(ts(pk,:,:),1));
         
h = mean(h_arr,2);
eh = std(h_arr,0,2)/sqrt(nsub);

h(1) = mean(h_arr(1,coh4),2);
eh(1) = std(h_arr(1,coh4),0,2)/sqrt(length(find(coh4)));

     
function data = ts_detrend(data,tps,order)

        sz = size(data);
        
            for j=1:size(data,2)
                for k=1:size(data,3)
                    ts = squeeze(data(tps,j,k));
                    p = polyfit(tps,ts',order);
                    t = 1:sz(1);
                    
                    bl = zeros(1,sz(1));
                    for ip=0:order
                       bl = bl+t.^(order-ip)*p(ip+1);
                    end
                    data(:,j,k) = data(:,j,k) - bl'; 
                end
            end
        
        
     
        
        
        