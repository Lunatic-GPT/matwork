function plot_single
global pk;
pk = 8:12;

dlist = {'/media/SS_FMRI1/MRI_Analysis/Peng/RDK_12_16_2008/afni/results', ...
         '/media/SS_FMRI1/MRI_Analysis/Hua/rdk_12_04_08/afni/results', ...
         '/media/SS_FMRI1/MRI_Analysis/Zong/rdk_12_01_08/afni/results', ...
         '/media/SS_FMRI1/MRI_Analysis/KELIRIS004_GeorgeSister/rdk/afni/results', ...
         };

     dlist = {'results_Peng','results_Zong','results_Hua','results_K004'};
     
[a,astd,h_V3a]=get_ts(dlist{4},'ts_MOG.mat');
plot_ts(a,astd,'MOG_S4');
plot_hght(h_V3a,zeros(size(h_V3a)),'MOG_S4');

[a,astd,h_V3a]=get_ts(dlist{4},'ts_rSP.mat');
plot_ts(a,astd,'rSP');
plot_hght(h_V3a,zeros(size(h_V3a)),'rSP');


function plot_hght(h,eh,mask)

yrange =[-1.2,1];
figure;
    if length(h) == 4
       errorbar([0.125,0.25,0.5,1],h,eh,'ko-','MarkerSize',6);
    else 
       errorbar([0.25,0.5,1],h,eh,'ko-','MarkerSize',6);
    end
set(gcf,'Units','inch','Position',[3,4,3.4,2.2]);
set(gca,'FontSize',14);
xlim([0,1.1]);
ylim(yrange);
xlabel('Coherence');
ylabel('Signal Change (%)');
set(gca,'TickLength',[0.03,0.03])
set(gca,'XTick',0:0.2:1);
set(gca,'YTick',-1:0.5:2);


%text(5,0,mask,'FontSize',14);
saveas(gcf,['hght_',mask,'.fig']);

function plot_ts(a,a_std,mask)
yrange =[-0.5,1.5];
figure;
symb = {'r','g','b','k'};
for i=1:size(a,2)
    if isempty(a_std)
         plot(1:24,a(:,i),symb{i});
    else 
        errorbar(1:24,a(:,i),a_std(:,i),symb{i});
    end
    hold on;
end
set(gcf,'Units','inch','Position',[3,4,3.4,2.2]);
set(gca,'FontSize',14);
rectangle('Position',[4.5,yrange(1),5,0.2],'FaceColor',[0.5,0.5,0.5]);
xlim([0,24.5]);
ylim(yrange);
xlabel('t (TR)');
ylabel('Signal Change (%)');
set(gca,'TickLength',[0.03,0.03]);
set(gca,'XTick',0:5:25);
set(gca,'YTick',-1:0.5:2);
text(5,0,mask,'FontSize',14);
saveas(gcf,['ts_',mask,'.fig']);


function [a,a_std,h] = get_ts(dlist,fname)
  
    load(fullfile(dlist,fname),'ts_mean','ts_std');
    
 a = eval('ts_mean');
 a_std =eval('ts_std');
 
 a = ts_detrend(a,[1:4,21:24],1);
 global pk;
h = mean(a(pk,:),1);



function data = ts_detrend(data,tps,order)

        sz = size(data);
        
            for j=1:size(data,2)
                    ts = squeeze(data(tps,j));
                    p = polyfit(tps,ts',order);
                    t = 1:sz(1);
                    
                    bl = zeros(1,sz(1));
                    for ip=0:order
                       bl = bl+t.^(order-ip)*p(ip+1);
                    end
                    data(:,j) = data(:,j) - bl'; 
            end
        
        
     
        
        
        