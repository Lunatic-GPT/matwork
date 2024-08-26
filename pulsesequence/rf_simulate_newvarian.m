function [mx,my,mz]=rf_simulate_newvarian(fname,pw,B1max,freq,M0)
% rf_simulate_varian(fname,pw,B1max,freq,M0)
% fname: name of the rf pattern
% pw : total pulse width in seconds.
% B1max: peak RF field strength in Hz. (Larmor frequency of proton at that field) 
% freq: frequencies to simulate. in Hz.
% M0: a 1*3 matrix specifying the initial magnetization. default:[0,0,1]



T1 = 2.01;
T2 = 0.2;

if length(B1max)==1
    B1max=B1max*ones(size(pw));
end


B1max = B1max/4258;  %convert to Gauss.
fname = str2cell(fname);

a=[];

if ~strcmp(computer,'PCWIN64')
    if exist('/home/xiaopeng/vnmrsys','dir')
      root='/home/xiaopeng/vnmrsys/shapelib';
    elseif exist('/home/xiaopeng/labhome/vnmrsys','dir');
      root='/home/xiaopeng/labhome/vnmrsys/shapelib';
    else
        root='/home/xiaopeng/vesta/vnmrsys/shapelib';
    end
     
else
    root = 'c:/labhome/vnmrsys/shapelib';
end

for i=1:length(pw)

    if isempty(fname{i})
        a_tmp=[0,0,1];
    elseif strcmp(fname{i},'x')
        a_tmp = [0,B1max(i),1];
    elseif strcmp(fname{i},'y')
        a_tmp = [pi/2,B1max(i),1];
    elseif strcmp(fname{i},'-x')
        a_tmp = [pi,B1max(i),1];
    elseif strcmp(fname{i},'-y')
        a_tmp = [pi*3/2,B1max(i),1];
    else
        a_tmp=textread(fullfile(root,fname{i}),'','commentstyle','shell');
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
   [mx,my,mz]=Bloch(rf,zeros(size(rf)),t,T1,T2,freq(:),0,0);
else
    M0 = repmat(M0(:),[1,length(freq)]);
   [mx,my,mz]=Bloch(rf,zeros(size(rf)),t,T1,T2,freq(:),0,0,M0(1,:),M0(2,:),M0(3,:));
end    
figure;subplot(2,1,1);
plot(freq,mz(:),'k','LineWidth',2);
hold on;plot(freq,sqrt( my(:).^2+mx(:).^2),'r','LineWidth',2);
legend('mz','mxy');
xlabel('Frequency (Hz)');

ph = atan2(my,mx)*180/pi;
subplot(2,1,2);plot(freq,ph,'g','LineWidth',1);
xlabel('Frequency (Hz)');
ylabel('phase');
ylim([-180,180]);


