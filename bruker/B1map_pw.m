function B1map_pw(f_prefix,pw90)
% change the pwr while fixing the pulse width
% number of peaks in the abs(singal) vs power plot.
%pw90 in units of us.
% B1map_pw(f_prefix,pw90)
 
if exist([f_prefix,'_recon+orig.HEAD'],'file')
  [d,info]=BrikLoad([f_prefix,'_recon+orig']);
else
  [d,info]=BrikLoad([f_prefix,'+orig']);
end

if ~exist('pw90','var')
pw90=1000;
end

refatt=readbPar([f_prefix,'/method'],'SLRefAtt');
if ndims(d) ==3
    d = reshape(d,[size(d,1),size(d,2),1,size(d,3)]);
end

%tpwr1=readbPar([f_prefix,'/method'],'PVM_ppgPowerList1');

%pw = readbPar([f_prefix,'/acqp'],'P');

%pw=pw(15)/1e6;
pw=readbPar([f_prefix,'/method'],'tSL');
pw=pw*1000000;
b1= zeros(size(d,1),size(d,2),size(d,3),2);
x=pw;
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

            f0=pi/2/pw90;
            [b,r] = nlinfit(x(:),y(:),@cosfita,[max(y),f0]);
            b1(i,j,k,1) = b(2)/2/pi*1e6; %convert to hertz at the maximum power used.
            b1(i,j,k,2)=1-sum(r.^2)/sum(y.^2);
            
        end
    end
end

b1=refatt+factor2db(b1(:,:,:,1)/250);
WriteBrikEZ(b1,info,'epi_B1map',[f_prefix,'_RefAtt_Map']);




