function batch_plot_single

sid = 'Hua'; 
view ='+orig';
add_rdk_subject(sid);
 %thr = 1.92;nvox = 31;
 %thr = 3.29;nvox = 7;
 %glm_mask(sid,thr,nvox);
 
 
% mstr = visual_areas(sid);
 %mstr = visual_areas_ecc;
%mstr = {'mask_MTLoc+orig;step(a)','V5u'};
%mstr = {'V5u_ap+tlrc;equals(a,1)','V5u_a';'V5u_ap+tlrc;equals(a,2)','V5u_p'; ...
%        'V5m_ap+tlrc;equals(a,1)','V5m_a';'V5m_ap+tlrc;equals(a,2)','V5m_p'};
mstr = glm_areas(sid);

if ~strcmp(sid,'K004')
flist = {['ts_sort_coh1',view],['ts_sort_coh2',view],['ts_sort_coh3',view],['ts_sort_coh4',view]};
else
 flist = {['ts_sort_coh1',view],['ts_sort_coh2',view],['ts_sort_coh3',view]};
end

exc_list={[],[],[],[]};
for i=1:size(mstr,1)
  tsroi_avetrl(flist,mstr{i,1},mstr{i,2},20,exc_list);
  plot_ts(mstr{i,2});
end

% note: Hua: some positive activation in V4 at 100%. 
%       Peng: small negative(positive) activation in V1(V2).positive
%       activation in V3 and V4.
function m=glm_areas(sid)

if strcmp(sid,'Zong')
    m(1,:) = {'mask_REML_coh234_thr3.29_clst+orig;or(equals(a,1),equals(a,3),equals(a,9))','V5glm'};
    m(2,:) = {'mask_REML_coh234_thr3.29_clst+orig;or(equals(a,2),equals(a,6))','STglm'};
    m(3,:) = {'mask_REML_coh234_thr3.29_clst+orig;or(equals(a,5),equals(a,8))','V3Aglm'};
    m(4,:) = {'mask_REML_coh234_thr3.29_clst+orig;or(equals(a,4),equals(a,11))','IPSglm'};
    m(5,:) = {'mask_REML_coh234_thr3.29_clst+orig;equals(a,10)','CCglm'};  %cingulate cortex
    m(6,:) = {'mask_REML_coh234_thr3.29_clst+orig;equals(a,7)','LSglm'};
elseif strcmp(sid,'Peng')
    m(1,:) = {'mask_REML_coh234_thr3.29_clst+orig;or(equals(a,1),equals(a,3))','V5glm'};
    m(2,:) = {'mask_REML_coh234_thr3.29_clst+orig;or(equals(a,2),equals(a,5),equals(a,8),equals(a,9))','V3Aglm'};
    m(3,:) = {'mask_REML_coh234_thr3.29_clst+orig;or(equals(a,4),equals(a,6),equals(a,7),equals(a,10),equals(a,14))','IPSglm'};
    m(4,:) = {'mask_REML_coh234_thr3.29_clst+orig;or(equals(a,11),equals(a,13))','LSglm'};
    m(5,:) = {'mask_REML_coh234_thr3.29_clst+orig;equals(a,12)','PreCglm'};
elseif strcmp(sid,'Hua')
    m(1,:) = {'mask_REML_coh234_thr3.29_clst+orig;or(equals(a,1),equals(a,2))','V5glm'};
    m(2,:) = {'mask_REML_coh234_thr3.29_clst+orig;equals(a,3)','STglm'};
    m(3,:) = {'mask_REML_coh234_thr3.29_clst+orig;equals(a,4)','V3Aglm'};
elseif strcmp(sid,'K004')
    m(1,:) = {'mask_REML_coh234_thr3.29_clst+orig;or(equals(a,1),equals(a,2))','V5glm'};
    m(2,:) = {'mask_REML_coh234_thr3.29_clst+orig;equals(a,3)','LSglm'};
    m(3,:) = {'mask_REML_coh234_thr3.29_clst+orig;equals(a,4)','STglm'};
    m(4,:) = {'mask_REML_coh234_thr3.29_clst+orig;equals(a,5)','V3Aglm'};
end
function plot_ts(mask)

load(['ts_',mask,'.mat'],'ts_mn','ts_sem');
%yrange =[-0.8,1.5];
if ~exist('ts_mn','var')
    return;
end
ts_mn = ts_mn*100;
ts_sem = ts_sem*100;
figure;
symb = {'r','g','b','k'};
for i=1:size(ts_mn,2)
    errorbar(1:size(ts_mn,1),ts_mn(:,i),ts_sem(:,i),symb{i});
    hold on;
end

set(gcf,'Units','inch','Position',[3,4,3.4,2.2]);
set(gca,'FontSize',14);
%rectangle('Position',[4.5,yrange(1),5,0.2],'FaceColor',[0.5,0.5,0.5]);
xlim([0,size(ts_mn,1)]);
yrange = min_max(cat(1,ts_mn+ts_sem,ts_mn-ts_sem));
ylim(yrange);
xlabel('t (TR)');
ylabel('Signal Change (%)');
set(gca,'TickLength',[0.03,0.03]);
set(gca,'XTick',0:5:25);
set(gca,'YTick',-1:0.5:2);
title(mask);
saveas(gcf,['ts_',mask,'.fig']);

  
%(uncorrected p=0.055, cluster size = 31 will give corrected p = 0.05);
%  mask_thr = 'mask_REML_coh3+orig';
  

function mstr = visual_areas(sid)

r = 'sroi2v_rh+orig';
l = 'sroi2v_lh+orig';
    
    
if strcmp(sid,'Zong') || strcmp(sid,'Hua') || strcmp(sid,'Peng')
    mstr(1,:) = {sprintf('%s[3];%s[2];%s[3];%s[2];or(a*equals(b,1),c*equals(d,1))',l,l,r,r),'V1u'};
    mstr(2,:) = {sprintf('%s[4];%s[5];%s[2];%s[4];%s[5];%s[2];or((a+b)*equals(c,1),(d+e)*equals(f,1))',l,l,l,r,r,r),'V2u'};
    mstr(3,:) = {sprintf('%s[6];%s[7];%s[2];%s[6];%s[7];%s[2];or((a+b)*equals(c,1),(d+e)*equals(f,1))',l,l,l,r,r,r),'V3u'};
    mstr(4,:) = {sprintf('%s[8];%s[2];%s[8];%s[2];or(a*equals(b,1),c*equals(d,1))',l,l,r,r),'V4u'};
    mstr(5,:) = {sprintf('%s[9];%s[2];%s[9];%s[2];or(a*equals(b,1),c*equals(d,1))',l,l,r,r),'V5u'};
    
    mstr(6,:) = {sprintf('%s[19];%s[19];or(a,b)',l,r),'V1m'};
    mstr(7,:) = {sprintf('%s[20];%s[21];%s[20];%s[21];or(a,b,c,d)',l,l,r,r),'V2m'};
    mstr(8,:) = {sprintf('%s[22];%s[23];%s[22];%s[23];or(a,b,c,d)',l,l,r,r),'V3m'};
    mstr(9,:) = {sprintf('%s[24];%s[24];or(a,b)',l,r),'V4m'};
    mstr(10,:) = {sprintf('%s[25];%s[25];or(a,b)',l,r),'V5m'};
    
    if strcmp(sid,'Hua')
      mstr(3,:) = {sprintf('%s[7];%s[2];%s[6];%s[7];%s[2];or(a*equals(b,1),(c+d)*equals(e,1))',l,l,r,r,r),'V3u'};
      mstr(8,:) = {sprintf('%s[23];%s[22];%s[23];or(a,b,c)',l,r,r),'V3m'};
    end
    
    if strcmp(sid,'Peng')
      mstr(3,:) = {sprintf('%s[7];%s[2];%s[6];%s[7];%s[2];or(a*equals(b,1),(c+d)*equals(e,1))',l,l,r,r,r),'V3u'};
      mstr(8,:) = {sprintf('%s[23];%s[22];%s[23];or(a,b,c)',l,r,r),'V3m'};
      mstr(2,:) = {sprintf('%s[5];%s[2];%s[4];%s[5];%s[2];or(a*equals(b,1),(c+d)*equals(e,1))',l,l,r,r,r),'V2u'};
      mstr(7,:) = {sprintf('%s[21];%s[20];%s[21];or(a,b,c)',l,r,r),'V2m'};
    end
else
    mstr(1,:) = {sprintf('%s[3];%s[2];%s[3];%s[2];or(a*equals(b,1),c*equals(d,1))',l,l,r,r),'V1u'};
    mstr(2,:) = {sprintf('%s[4];%s[5];%s[2];%s[4];%s[5];%s[2];or((a+b)*equals(c,1),(d+e)*equals(f,1))',l,l,l,r,r,r),'V2u'};
    mstr(3,:) = {sprintf('%s[6];%s[7];%s[2];%s[6];%s[7];%s[2];or((a+b)*equals(c,1),(d+e)*equals(f,1))',l,l,l,r,r,r),'V3u'};
    mstr(4,:) = {sprintf('%s[8];%s[2];%s[8];%s[2];or(a*equals(b,1),c*equals(d,1))',l,l,r,r),'V4u'};
    mstr(5,:) = {sprintf('%s[17];%s[17];or(a,b)',l,r),'V1m'};
    mstr(6,:) = {sprintf('%s[18];%s[19];%s[18];%s[19];or(a,b,c,d)',l,l,r,r),'V2m'};
    mstr(7,:) = {sprintf('%s[20];%s[21];%s[20];%s[21];or(a,b,c,d)',l,l,r,r),'V3m'};
    mstr(8,:) = {sprintf('%s[22];%s[22];or(a,b)',l,r),'V4m'};
end

function mstr = visual_areas_ecc


l = 'roi_eccgt180_lh+orig';
r = 'roi_eccgt180_rh+orig';
mstr = cell(8,2);
  for i=1:4
   % mstr{(i-1)*2+1,1} = sprintf('V%du+orig;%s[3];%s[2];%s[3];%s[2];and(a,b*equals(c,1)+d*equals(e,1))',i,l,l,r,r);
   % mstr{(i-1)*2+1,2} = sprintf('V%du_le',i);
   % mstr{(i-1)*2+2,1} = sprintf('V%du+orig;%s[4];%s[2];%s[4];%s[2];and(a,b*equals(c,1)+d*equals(e,1))',i,l,l,r,r);
   % mstr{(i-1)*2+2,2} = sprintf('V%du_se',i);
    mstr{(i-1)*2+1,1} = sprintf('V%dm+orig;%s[3];%s[2];%s[3];%s[2];and(a,b*equals(c,1)+d*equals(e,1))',i,l,l,r,r);
    mstr{(i-1)*2+1,2} = sprintf('V%dm_le',i);
    mstr{(i-1)*2+2,1} = sprintf('V%dm+orig;%s[4];%s[2];%s[4];%s[2];and(a,b*equals(c,1)+d*equals(e,1))',i,l,l,r,r);
    mstr{(i-1)*2+2,2} = sprintf('V%dm_se',i);
  end
    
    
    
    
    
    
 
 
