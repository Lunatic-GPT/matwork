function epi_B1map_pwr(f_prefix)
% change the pwr while fixing the pulse width
[d,info]=BrikLoad([f_prefix,'+orig']);
if ndims(d) ==3
    d = reshape(d,[size(d,1),size(d,2),1,size(d,3)]);
end

tpwr1=readPar([f_prefix,'.fid/procpar'],'tpwr1');
ref_pos=readPar([f_prefix,'.fid/procpar'],'ref_pos');
pw = readPar([f_prefix,'.fid/procpar'],'p1');

tpwr1(ref_pos+1)=[];

[t_sort,ind]=sort(tpwr1);
   
d_sort = d(:,:,:,ind);

WriteBrikEZ(d_sort,info,'epi_B1map',[f_prefix,'_sort']);

b1= zeros(size(d,1),size(d,2),size(d,3),2);
x=10.^(tpwr1/20);
            
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
            [b,r] = nlinfit(x(:),y(:),@sinfita,[max(y),f0]);
            b1(i,j,k,1) = b(2)*x(end)/pw/2/pi; %convert to hertz at the maximum power used.
            b1(i,j,k,2)=1-sum(r.^2)/sum(y.^2);
        end
    end
end

WriteBrikEZ(b1,info,'epi_B1map',[f_prefix,'_map']);

