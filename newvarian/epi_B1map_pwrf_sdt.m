function epi_B1map_pwrf_sdt(f_prefix)
% epi_B1map_pwrf(f_prefix)
% change the fine pwr while fixing the pulse width
d=rdSdt(f_prefix);

tpwr1=readPar([f_prefix,'.fid'],'tpwr1f');
img=readPar(f_prefix,'image');
tpwr1=tpwr1(img==1);

pw = readPar([f_prefix,'.fid'],'p1');

[t_sort,ind]=sort(tpwr1);
   
d_sort = d(:,:,:,ind);

writesdt4(d_sort,[f_prefix,'_sort']);
b1= zeros(size(d,1),size(d,2),size(d,3),2);

x = t_sort;
d=d_sort;

dmean = mean(d,4);
dmax = max(dmean(:));
rfpat = readPar([f_prefix,'.fid'],'p1pat');

if strcmp(rfpat,'"sinc"')
    factor = 0.1777;
elseif strcmp(rfpat,'"SGLsinc"');
    factor = 0.1777;
else
    error('rf pattern unknown');
end


options=optimset('MaxIter',50,'Display','off');
for i=1:size(d,1)
    fprintf('%d  ',i);
    for j=1:size(d,2)
        for k=1:size(d,3)
            y=d(i,j,k,:);
            y=y(:);
            
            if max(y)<dmax/20
                continue;
            end
            
            maxtab = peakdet(y,max(y)/50);
            if isempty(maxtab)
                k0 = 2*pi/max(x)/4;
            else
                k0 = 2*pi/x(maxtab(1,1))/5;
            end
            [b,r] = lsqcurvefit(@sinfita,[max(y),k0],x(:),y(:),[],[],options);
            if i==20 && j==22
                disp('pause');
            end
            b1(i,j,k,1) = b(2)*4095/pw/2/pi/factor; %convert to hertz at the maximum power used for a hard pulse.
            b1(i,j,k,2)=1-sum(r.^2)/sum(y.^2);
        end
    end
end
writesdt4(b1(:,:,:,1),[f_prefix,'_B1map']);

