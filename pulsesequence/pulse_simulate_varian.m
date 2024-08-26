function [mx,my,mz]=pulse_simulate_varian(rfname,pw,B1max,linewidth,freq0,M0)
% [mx,my,mz]=pulse_simulate_varian(rfname,pw,B1max,linewidth,freq0,M0)
% rfname: name of the rf pattern.
% pw : total pulse width in seconds.
% B1max: peak RF field strength in Hz. (Larmor frequency of proton at that field)
% linewidth: FWHM of a Gaussian shaped spectrum. frequencies to simulate is -4*linewidth to 4*linewidth. in Hz. 
% freq0: center resonance frequency.  (default zero)
% M0: a 1*3 matrix specifying the initial magnetization. default:[0,0,1]


T1 = 1.01;
T2 = 0.2;
rfname = str2cell(rfname);

if ~strcmp(computer,'PCWIN64')
    if exist('/data/users/xiaopeng/vnmrsys','dir')
      root='/data/users/xiaopeng/vnmrsys/shapelib';
    elseif exist('/home/xiaopeng/labhome/vnmrsys','dir');
      root='/home/xiaopeng/labhome/vnmrsys/shapelib';
    else
        root='/home/xiaopeng/vesta/vnmrsys/shapelib';
    end
     
else
    root = 'z:/home/xiaopeng/vnmrsys/shapelib';
end

if ~exist('freq0','var')
    freq0= 0;
end

if linewidth >0
    freq = linspace(-linewidth*4,linewidth*4,800)+freq0;
else
    freq = freq0;
end

if length(B1max)==1
    B1max=B1max*ones(length(rfname));
end

if length(pw)==1
    pw=pw*ones(length(rfname));
end

B1max = B1max/4258;  %convert to Gauss.

a=[];
for i=1:length(pw)

    if isempty(rfname{i})
        a_tmp=repmat([0,0,1],[100,1]);     
    elseif strcmp(rfname{i},'x')
        a_tmp = [-pi/2,B1max(i),1];
    elseif strcmp(rfname{i},'y')
        a_tmp = [0,B1max(i),1];
    elseif strcmp(rfname{i},'-x')
        a_tmp = [pi/2,B1max(i),1];
    elseif strcmp(rfname{i},'-y')
        a_tmp = [pi,B1max(i),1];
    else
        a_tmp=textread(fullfile(root,rfname{i}),'','commentstyle','shell');
        a_tmp(:,1)=a_tmp(:,1)*pi/180;
        a_tmp(:,2) = a_tmp(:,2)/max(a_tmp(:,2))*B1max(i);
    end
    
 if i==1
     a=a_tmp(:,1:2);
     t = a_tmp(:,3);
     t = t/sum(t)*pw(i);
 else
     a=cat(1,a,a_tmp(:,1:2));
     t_tmp = a_tmp(:,3);
     t = cat(1,t,t_tmp/sum(t_tmp)*pw(i));
 end

end
 
rf = a(:,2).*(cos(a(:,1))+1i*sin(a(:,1)));

if ~exist('M0','var')
   [mx,my,mz]=bloch(rf,zeros(size(rf)),t,T1,T2,freq(:),0,2);
else
    M0 = repmat(M0(:),[1,length(freq)]);
   [mx,my,mz]=bloch(rf,zeros(size(rf)),t,T1,T2,freq(:),0,2,M0(1,:),M0(2,:),M0(3,:));
end    


if linewidth >0
    
figure;subplot(2,1,1);
plot(freq,mz(end,:),'k','LineWidth',2);
hold on;plot(freq,sqrt( my(end,:).^2+mx(end,:).^2),'r','LineWidth',2);
legend('mz','mxy');
xlabel('Frequency (Hz)');

ph = atan2(my,mx)*180/pi;
subplot(2,1,2);plot(freq,ph(end,:),'g','LineWidth',1);
xlabel('Frequency (Hz)');
ylabel('phase (degree)');
ylim([-180,180]);

figure;
sig = linewidth/2/1.17741;
a=gauss(0,sig,freq);
a=repmat(a,[size(mx,1),1])/sum(a);
mx = mx.*a;
my=my.*a;
mz=mz.*a;


mxy_total = sqrt(sum(mx,2).^2+sum(my,2).^2);
mz_total = sum(mz.*a,2);

plot(cumsum(t),mz_total,'k',cumsum(t),mxy_total,'b');
xlabel('Time (s)');
ylabel('Mxy/Mz');
legend('Mz','Mxy');

end