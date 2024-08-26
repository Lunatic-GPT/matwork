function [ph,trig]=physioTriggerNew(data,minInterval,t,nphase,trig_delay,thr,plot_results,fakePeaks,Peaks2MON,realPeaks,interpPeaks)
% [ph,trig]=physioTrigger(data,minInterval,TR,plot_results)
% data: physiology signal
% minInterval: unit: s; the minimum spacing between peaks
% TR: unit: s
% nphase: number of phases to divide each cardiac cycle into
% trig_delay: number of time points to shift the trig for defining the
% phase intervals.  Can be non-integer.

%output:
%ph: phase indices for each physio signal time point: 1 to nphase
% trig: peak positions for the physio signal

if ~exist('plot_results','var')
    plot_results=false;
end

if ~exist('realPeaks','var')
    realPeaks= [];
end

fakePeaks=cat(1,fakePeaks(:),Peaks2MON(:));

data=double(data);
%nTRmin=round(minInterval/TR);

mx_orig=peakdet(data,100);
i_orig=mx_orig(:,1);
%% add real peaks;
i_real=ind_match(t,realPeaks);
i_orig=sort([i_orig(:);i_real(:)]);

%% the first point is a fake peak
if i_orig(1)==1
    i_orig(1)=[];
end

%% remove fake peaks

%[mx,imx]=setdiff(roundp1(t(i_orig)),roundp1(fakePeaks));

ii_rm=ind_match(t(i_orig),fakePeaks);
i_new=i_orig;
i_new(ii_rm)=[];


% 
% if ~isempty(fakePeaks)
%     tround=roundp1(t(i_orig));
%     ffakePeaks=setdiff(roundp1(fakePeaks),tround);
%     
%     if length(ffakePeaks)>0
%         fprintf('Some fake peaks not detected\n');
%         disp(ffakePeaks);
%     end
% end


%% remove peaks too close to real peaks.

%% ends here
%if isempty(ind_rm)
 %   break;
%end
%end

interval2=diff(t(mx(:,1)));

interval2_old=interval2;
sd=std(interval2);

mn=mean(interval2);
if thr>0
thr=thr*sd;
ind_outlier=find(abs(interval2-mn)>thr);
i=1;

ind_rm=[];
while i<length(ind_outlier)
    
    if ind_outlier(i+1)==ind_outlier(i)+1  % there is a fake peak or noisy peak
        if sign(interval2(ind_outlier(i))-mn)*sign(interval2(ind_outlier(i+1))-mn)==-1  % one peak position was off
            interval_tmp=interval2(ind_outlier(i))+interval2(ind_outlier(i+1));
            interval2(ind_outlier(i))=round(0.5*interval_tmp);
            interval2(ind_outlier(i)+1)=interval_tmp-interval2(ind_outlier(i));  
           
        elseif sign(interval2(ind_outlier(i))-mn)==-1 && sign(interval2(ind_outlier(i+1))-mn)==-1% combine two peaks into one; one peak was fake
            
            interval2(ind_outlier(i+1))=interval2(ind_outlier(i))+interval2(ind_outlier(i+1));
            
            ind_rm(end+1)=ind_outlier(i);
        else
            warning('Not sure what to do');
        end    
    
         i=i+2;
    elseif interval2(ind_outlier(i))>mn % there is arrhythmia
        
        i=i+1;
    else
        warning('This should not happen');
    end
end

else
   ind_rm=[]; 
end

interval2(ind_rm)=[];

mx2=mx;
mx2(ind_rm,:)=[];
%tmp=cumsum(interval2)+t(mx2(1,1));

tmx2=[t(mx2(:,1))',mx2(:,2)];

if ~isempty(Peaks2MON)
    
    for i=1:length(Peaks2MON)
[tmp,ind]=min(abs(t-roundp1(Peaks2MON(i))));
tmx2=cat(1,tmx2,[roundp1(Peaks2MON(i))',data(ind)]);
mx2=cat(1,mx2,[ind,data(ind)]);
    end
end

tmx2(:,1)=roundp1(tmx2(:,1));

[tmp,ind]=sort(tmx2(:,1));
tmx2=tmx2(ind,:);
mx2=mx2(ind,:);


% tmpPeaks=setdiff(Peaks2MON,mx2(:,1));
% if length(tmpPeaks)>0
%     fprintf('Some Peaks 2 MON not detected\n');
%     disp(tmpPeaks);
% 
% end

tmx3=tmx2;

mx3=mx2;

for i=1:length(Peaks2MON)
[tmpPeaks,imx3]=intersect(tmx3(:,1),roundp1(Peaks2MON(i)));

tmp=((tmx3(imx3+1,1)+tmx3(imx3-1,1))/2);

[tmp,tind]=min(abs(tmp-t));

if ~isempty(tmp)
mx3(imx3,1)=tind;
mx3(imx3,2)=data(mx3(imx3,1));
end

end

trig=t(mx3(:,1));
ph = getIndexPhysioSignalPC(trig,t,length(data),nphase,trig_delay);

%fprintf('percent of interval less than phase number %d - %f\n',nphase,sum(diff(trig)<nphase*tr*2)/length(trig));


%%
if plot_results
    %%
   
    tt=max(t);
    nseg=ceil(max(t)/30);
    
    clr=jet(nphase);
    
    [peak_real,imx2]=intersect(mx2(:,1),mx_orig(:,1));
    peak_real(:,2)=mx2(imx2,2);
    
    [peak_real_adjust,imx3]=setdiff(mx3(:,1),mx2(:,1));
    peak_real_adjust(:,2)=mx3(imx3,2);
    
    [peak_fake,imx]=setdiff(mx_orig(:,1),mx2(:,1));
    peak_fake(:,2)=mx_orig(imx,2);
    
    for i=1:nseg
        
      h=figure(100+i);%   subplot(3,ceil(nseg/3),i);
        t_thr1=tt/nseg*(i-1);
        t_thr2=tt/nseg*i;
        
        
      
        ind=find(t<=t_thr2&t>t_thr1);
        
%         for j=1:length(ind)
%           plot(t(ind(j)),data(ind(j)),'.','Color',clr(ph(ind(j)),:));
%           hold on;
%         end
        
      plot(t(ind),data(ind),'k.');
      hold on;
      plot(t(ind),ph(ind)/nphase*diff(min_max(data))+double(max(data)));
      
        hold on;
        ind=t(peak_real(:,1))<=t_thr2&t(peak_real(:,1))>t_thr1;
         ind3=t(peak_real_adjust(:,1))<=t_thr2&t(peak_real_adjust(:,1))>t_thr1;
         
        ind2=t(peak_fake(:,1))<=t_thr2&t(peak_fake(:,1))>t_thr1;
        
      hp1=  plot(t(peak_real(ind,1)),peak_real(ind,2),'ob','MarkerSize',4);
       hp2= plot(t(peak_fake(ind2,1)),data(peak_fake(ind2,1)),'*r','MarkerSize',4);
       hp3= plot(t(peak_real_adjust(ind3,1)),peak_real_adjust(ind3,2),'or','MarkerSize',4);
       
      hp=hp1;
      lbl={'Real'};
      
      if ~isempty(hp2)
          hp(end+1)=hp2;
          lbl{end+1}='Fake';
      end
      
      if ~isempty(hp3)
          hp(end+1)=hp3;
          lbl{end+1}='MON';
      end
      
       legend(hp,lbl);
       
       
        xlabel('Time (s)');
        xlim([t_thr1,t_thr2]);
        
      set(h,'Position', [165         463        1475         420]);
      drawnow;
      hold off;
    end
 figure(200);
 
 interval3=diff(t(mx3(:,1)));
 
        plot(t(mx3(2:end,1)),interval3,'ro-');
        hold on;
   %     plot(mx(2:end,1)*TR,interval2_old*TR,'bo-');
        
        plot([0,max(t)],[mn-thr,mn-thr],'k-');
        
        plot([0,max(t)],[mn+thr,mn+thr],'k-');
        
        xlabel('Time (TR)');
        xlim([0,max(t)]);
        ylabel('Interval (s)');
        
        title(sprintf('Heart Rate = %3.1f/min',60/mn));
        
        hold off;

end

function i_new = rm_peaks_too_close_to_real(t,i_orig,i_real,minInterval)


ind_rm=[];

for i=1:length(i_real)
   
    
interval=roundp1(t(i_))-roundp1(realPeaks(i));


ind_tmp=find(abs(interval)<minInterval & abs(interval)~=0);

ind_rm(end+1:end+length(ind_tmp))=ind_tmp;

end

mx(ind_rm,:)=[];


%while 1
interval=diff(t(mx(:,1)));

i2=find(interval<minInterval);

ind_rm=[];
j=1;


while j<=length(i2)
   
    if i2(j)==1
        
        if   mx(1,1)==1  %the first peak is a fake peak
            ind_rm(end+1)=1;
        end
        j=j+1;
    
    else
        
 %old method; works in HVB320 and HVB319; but not in HVB321      
%         if interval(i2(j)-1)<interval(i2(j)+1)  % the following interval is normal
%             ind_rm(end+1)=i2(j);
%         else         % the previous interval is normal; 
%                      %if a fake peak is present, in_rm generated by the two conditions may overlap
%             ind_rm(end+1)=i2(j)+1;
%         end
        
       if j<length(i2) && i2(j)+1 ==i2(j+1) % the next interval also too short, then combine the two; 
           
           ind_rm(end+1)=i2(j)+1;
           
           j=j+2;
       else  % only itself is too short
            ind_rm(end+1)=i2(j)+1;
           
           j=j+1;
       end
    end
    
end
mx(ind_rm,:)=[];


function res=ind_match(t,realPeak)

for i=1:length(realPeak)
    
    [tmp,res(i)]=min(abs(t-realPeak(i)));
    
end
    
    
function ph = getIndexPhysioSignalPC(trig,t,npts,nphase,trig_delay)
 
% trig: time for peaks;  
% npts: number of physio signal time points
% nphase: number of phases to divide each cardiac cycle into

ph=zeros(1,npts);
trig=trig+trig_delay;
for i=1:npts
 
    itr=find(t(i)>=trig(1:end-1)&t(i)<trig(2:end));
    
  
    if isempty(itr)
        if t(i)<trig(1) %before first peak
      
            trig0=2*trig(1)-trig(2);
            
          ph(i)=round((t(i)-trig0)/(trig(2)-trig(1))*nphase);
          
        else  %after last peak
        
           if round((t(npts)-trig(end))/(trig(end)-trig(end-1))*nphase)>nphase % too long
               ph(i)=round((t(i)-trig(end))/(t(npts)-trig(end))*nphase);
           else
             ph(i)=round((t(i)-trig(end))/(trig(end)-trig(end-1))*nphase);
           end
        end
    else
      
        ph(i)=round((t(i)-trig(itr))/(trig(itr+1)-trig(itr))*nphase);
        
    end
    
end

ph(ph==0)=nphase;
for i=1:nphase
    fprintf('Phase %d = %d\n',i,sum(ph==i));
end


function f=roundp1(f)

f=round(f*100)/100;
f=round(f*10)/10;




