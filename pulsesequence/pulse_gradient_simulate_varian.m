function [mx,my,mz]=pulse_gradient_simulate_varian(rfname,pw,B1max,gr,dvox,tof,M0)
% [mx,my,mz]=pulse_gradient_simulate_varian(rfname,pw,B1max,gradient,dvox,tof,M0)
% rfname: name of the rf pattern. M*1
% pw : total pulse width in seconds. 
% B1max: peak RF field strength in Hz. (Larmor frequency of proton at that field)
% gradient: gradient M*1, 2, or 3. in G/cm
% dvox: voxel size in cm.
% tof: resonance frequency offset, Unit: Hz.
% M0: a 1*3 matrix specifying the initial magnetization. default:[0,0,1]

if nargin == 0 
    help pulse_gradient_simulate_varian;
end

T1 = 1.01;
T2 = 0.2;
rfname = str2cell(rfname);

if size(gr,1)==1
    gr = gr';
end

steps = 10;  % number of points in each voxel to simulate.
nvox = 32;

dp = linspace(-(nvox-0.5)*dvox,(nvox+0.5)*dvox,nvox*steps);
  
if size(gr,2)==1
  dp = dp';
elseif size(gr,2)==2
  [dp1,dp2] = meshgrid(dp,dp);
  dp = [dp1(:),dp2(:)]*dvox;
elseif size(gr,2)==3
   [dp1,dp2,dp3] = meshgrid(dp,dp,dp);
   dp=[dp1(:),dp2(:),dp3(:)]*dvox;
end

if ~exist('tof','var')
    tof =0;
end

if length(B1max)==1
    B1max=B1max*ones(length(rfname));
end

if length(pw)==1
    pw=pw*ones(length(rfname));
end

B1max = B1max/4258;  %convert to Gauss.

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
        a_tmp=textread(rfname{i},'','commentstyle','shell');
        a_tmp(:,1)=a_tmp(:,1)*pi/180;
        a_tmp(:,2) = a_tmp(:,2)/max(a_tmp(:,2))*B1max(i);
    end
    
    gr2_tmp = repmat(gr(i,:),[size(a_tmp,1),1]);
        
 if i==1
     a=a_tmp(:,1:2);
     t = a_tmp(:,3);
     t = t/sum(t)*pw(i);
     gr2 = gr2_tmp;
 else
     a=cat(1,a,a_tmp(:,1:2));
     t_tmp = a_tmp(:,3);
     t = cat(1,t,t_tmp/sum(t_tmp)*pw(i));
     gr2 = cat(1,gr2,gr2_tmp);
 end

end
 
rf = a(:,2).*(cos(a(:,1))+1i*sin(a(:,1)));

if ~exist('M0','var')
   [mx,my,mz]=bloch(rf,gr2,t,T1,T2,tof,dp,0);
else
    M0 = repmat(M0(:)',[size(dp,1),1]);
   [mx,my,mz]=bloch(rf,gr2,t,T1,T2,tof,dp,0,M0(:,1),M0(:,2),M0(:,3));
end    

figure;subplot(2,1,1);
plot(dp(:,1),mz,'k','LineWidth',2);
hold on;plot(dp(:,1),sqrt( my.^2+mx.^2),'r','LineWidth',2);
legend('mz','mxy');
xlabel('Position (cm)');

ph = atan2(my,mx)*180/pi;
subplot(2,1,2);plot(dp(:,1),ph,'g','LineWidth',1);
xlabel('Position (cm)');
ylabel('phase');
ylim([-180,180]);



figure;

mz = reshape(mz,[steps,nvox]);
mxy = mx+1i*my;
mxy = reshape(mxy,[steps,nvox]);

mxy_total = mean(mxy,1);
mz_total = mean(mz,1);


dp_tmp = linspace(-(nvox-0.5)*dvox,(nvox+0.5)*dvox,nvox);

plot(dp_tmp,mz_total,'k',dp_tmp,abs(mxy_total),'b');
xlabel('Position (cm)');
ylabel('Mxy/Mz');
legend('Mz','Mxy');



