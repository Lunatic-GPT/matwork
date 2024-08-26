function epi_B1map_pwrf(f_prefix)
% epi_B1map_pwrf(f_prefix)
% change the fine pwr while fixing the pulse width
[d,info]=BrikLoad([f_prefix,'+orig']);
if ndims(d) ==3
    d = reshape(d,[size(d,1),size(d,2),1,size(d,3)]);
end

tpwr1=readPar([f_prefix,'.fid'],'tpwr1f');
ref_pos=readPar([f_prefix,'.fid'],'ref_pos');
pw = readPar([f_prefix,'.fid'],'p1');

tpwr1(ref_pos+1)=[];
d(:,:,:,ref_pos+1)=[];
[t_sort,ind]=sort(tpwr1);
   
d_sort = d(:,:,:,ind);

WriteBrikEZ(d_sort,info,'epi_B1map_pwrf',[f_prefix,'_sort']);

b1= zeros(size(d,1),size(d,2),size(d,3),2);

x = t_sort;
d=d_sort;

dmean = mean(d,4);
dmax = max(dmean(:));
rfpat = readPar([f_prefix,'.fid'],'p1pat');

if strcmp(rfpat,'"sinc"')
    factor = 0.1796;
elseif strcmp(rfpat,'"sinc5lapp"');
    factor = 0.1665;
else
    error('rf pattern unknown');
end


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
                k0 = 2*pi/max(x)/4;
            else
                k0 = 2*pi/x(maxtab(1,1))/4;
            end
            [b,r] = nlinfit(x(:),y(:),@sinfita,[max(y),k0]);
            if i==31 && j==37 && k==1
              disp('Pause');
            end
            b1(i,j,k,1) = b(2)*4095/pw/2/pi/factor; %convert to hertz at the maximum power used for a hard pulse.
            b1(i,j,k,2)=1-sum(r.^2)/sum(y.^2);
        end
    end
end

WriteBrikEZ(b1,info,'epi_B1map',[f_prefix,'_map']);

