function plot_motion_afni_mat_msub(dfile,ttl)
% dfile is a mat file;
% combine plots from different subjects into one figure for easy view
nsub=length(dfile);
[nrow,ncol]=nrow_ncol(nsub);
figure;

set(gcf,'Position',[1 750 ncol*225 nrow*240]);
for i=1:length(dfile)
    d=load(dfile{i});
    
    d=cat(2,d.an',squeeze(d.shift([3,1,2],1,:))');
    do_plot(d,nsub,i,ttl{i});
    
end


function do_plot(d,nsub,isub,ttl)

lbl={'IS','RL','AP','IS','RL','AP'};
% Test: the new image shifted to I, L, P
% Motion parameters were: S; R, and A; so the motion parameters are the amount of motion from new to base.
%
clr=lines(3);
[nrow,ncol]=nrow_ncol(nsub);

irow=ceil(isub*2/ncol);
icol=mod(isub*2-2,ncol)+1;

mg=0.02;
mgh=0.09;
hgt=(1-mgh)/nrow;
wdt=(1-mg)/ncol;
subplot('position',[mg+wdt*(icol-1),mgh+hgt*(nrow-irow),wdt-mg,hgt-mgh]);
for i=1:3
    plot(d(:,3+i),'.-','Color',clr(i,:));
    hold on;
end
set(gca,'FontSize',10);
%ylabel('(mm)');
xlim([0.5,size(d,1)+0.5]);
title([ttl,' - Translation']);
%legend(lbl(1:3));
set(gca,'xtick',0:20:160);
set(gca,'ytick',-3:0.2:3);

icol=mod(isub*2-1,ncol)+1;

subplot('position',[mg+wdt*(icol-1),mgh+hgt*(nrow-irow),wdt-mg,hgt-mgh]);

for i=1:3
    plot(d(:,i),'.-','Color',clr(i,:));
    hold on;
end
set(gca,'FontSize',10);
%ylabel('(degree)');
xlim([0.5,0.5+size(d,1)]);
title([ttl,' - Rotation']);
%legend(lbl(1:3));
%xlabel('TR index');
set(gca,'xtick',0:20:160);
set(gca,'ytick',-3:0.2:3);

set(gcf,'Units','pixels');


function [nrow,ncol]=nrow_ncol(nsub)
nrow=ceil(nsub/4);
ncol=nsub*2;
if nsub>4
    ncol=8;
end


%savetiffc(prefix);
