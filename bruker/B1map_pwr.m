function B1map_pwr(f_prefix,npeaks)
% change the pwr while fixing the pulse width
% number of peaks in the abs(singal) vs power plot.
[d,info]=BrikLoad([f_prefix,'+orig']);
if ndims(d) ==3
    d = reshape(d,[size(d,1),size(d,2),1,size(d,3)]);
end

tpwr1=readbPar([f_prefix,'/method'],'PVM_ppgPowerList1');

pw = readbPar([f_prefix,'/acqp'],'P');

pw=pw(15)/1e6;

b1= zeros(size(d,1),size(d,2),size(d,3),2);
x=10.^(-tpwr1/20);
disp(['Max power is ',num2str(min(tpwr1))]);            
[x,ind] = sort(x);
d=d(:,:,:,ind);   
dmax = max(d(:));
for i=1:size(d,1)
    for j=1:size(d,2)
        for k=1:size(d,3)
            y=d(i,j,k,:);
            y=y(:);
            
            if max(y)<dmax/20
                continue;
            end
            
            maxtab = peakdet(y,max(y)/50);
            if isempty(maxtab)
                f0 = 2*pi/max(x)/4;
            else
                f0 = 2*pi/x(maxtab(1,1))/4;
            end
            f0=2*pi/max(x)*npeaks/2;
            [b,r] = nlinfit(x(:),y(:),@sinfita,[max(y),f0]);
            b1(i,j,k,1) = b(2)*x(end)/pw/2/pi; %convert to hertz at the maximum power used.
            b1(i,j,k,2)=1-sum(r.^2)/sum(y.^2);
            
        end
    end
end

WriteBrikEZ(b1,info,'epi_B1map',[f_prefix,'_map']);

