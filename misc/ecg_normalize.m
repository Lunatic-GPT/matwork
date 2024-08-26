function ecg_normalize(fpattern,r_manual,vmin)
% ecg_normalize(fpattern[,r_manual,vmin])
% r_manual:  minimum value change to qualify as a new peak.
% if not set, will use 1/5 of the full data range for each file.
% vmin: minimum value to qualify as a peak.

dir_str = dir(fpattern);

for i=1:length(dir_str)
    a = load(dir_str(i).name);
    
    if exist('r_manual','var')
        r = r_manual;
    else
        r = (max(a)-min(a))/5;
    end
    
    [maxtab,mintab] = peakdet(a,r);
    
    if ~exist('vmin','var')
       vmin = (max(a)-min(a))/3+min(a);    
    end
      ind = maxtab(:,2)<vmin;
      maxtab(ind,:) = [];
      %mintab(ind,:) = [];
    
    
    fprintf('%s: hear rate %3.2f/min\n',dir_str(i).name,size(maxtab,1)*60*40/length(a));
    
    b = zeros(size(a));
    
    dip_r = find(mintab(:,1)>maxtab(1,1),1,'first');
    b(1:mintab(dip_r,1)) = a(1:mintab(dip_r,1))*1000/maxtab(1,2);
    dip_l = find(mintab(:,1)<maxtab(end,1),1,'last');
    b(mintab(dip_l,1):end) = a(mintab(dip_l,1):end)*1000/maxtab(end,2);
    
    for j=2:size(maxtab,1)-1
        pos = maxtab(j,1);
        dip_r = find(mintab(:,1)>pos,1,'first');
        dip_l = find(mintab(:,1)<pos,1,'last');
        b(mintab(dip_l):mintab(dip_r)) = a(mintab(dip_l):mintab(dip_r))*1000/maxtab(j,2);
            
    end
    
    
prefix = strtok(dir_str(i).name,'.');
save_mat_int(int16(b),[prefix,'_norm.1D']);


figure;
subplot(2,1,1);
plot(a);

t = strrep(dir_str(i).name,'_','\_');
title(t);
xlim([0,length(a)]);
subplot(2,1,2);
plot(b);
xlim([0,length(a)]);

set(gcf,'Units','pixels','Position',[6,653,1590,461]);
end
 
