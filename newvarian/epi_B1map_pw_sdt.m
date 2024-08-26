function epi_B1map_pw_sdt(f_prefix)
% epi_B1map_pwrf(f_prefix)
% change the pulse width pmt while fixing power
d=rdSdt(f_prefix);

img=readPar(f_prefix,'image');

pw = readPar([f_prefix,'.fid'],'pwc');
pw=pw(img==1);

[t_sort,ind]=sort(pw);
   
d_sort = d(:,:,:,ind);

writesdt4(d_sort,[f_prefix,'_sort']);
b1= zeros(size(d,1),size(d,2),size(d,3));

x = t_sort;
d=d_sort;

dmean = mean(d,4);
dmax = max(dmean(:));


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
            
            [~,ind_min]=min(y);
          
           % k0 = 1/(x(ind_min)*4)*2*pi;
           k0=3;
            b = lsqcurvefit(@cosfita,[max(y),k0],x(:),y(:),[],[],options);
            
            b1(i,j,k,1) = b(2)/2/pi; %convert to hertz at the maximum power used for a hard pulse.
            
        end
    end
end
writesdt4(b1(:,:,:,1),[f_prefix,'_B1map']);

