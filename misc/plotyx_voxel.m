function plotyx_voxel(fid_prefix,ijk,xpar,ind_ex,fit_func,beta0)
%plotyx_voxel(fid_prefix,ijk,xpar,ind_ex[,fit_func,beta0])
i=ijk(1)+1;
j=ijk(2)+1;
k=ijk(3)+1;


d= BrikLoadf([fid_prefix,'+orig']);

x=readPar(fid_prefix,xpar);
if size(d,4)>length(x)
  x = parValArray(fid_prefix,xpar);
end

ts = d(i,j,k,:);
ind = 1:size(d,4);
ind=setdiff(ind,ind_ex);

ts = ts(ind);
x=x(ind);

figure;
w = size(ts,1);
h = size(ts,2);
h2 = size(ts,3);
for i1=1:w
    for i2=1:h
        for i3=1:h2
            
            
          ind = ((i3-1)*h+i2-1)*w+i1;
          
          subplot(h*h2,w,ind);
          y = squeeze(ts(i1,i2,:));
          x = x(:);
          
          plot(x,real(y),'ko');
          hold on;          
          if any(~isreal(y))
           plot(x,imag(y),'k+');
          end
          ttl = sprintf('(%d,%d,%d)',i(i1)-1,j(i2)-1,k(i3)-1);
          title(ttl);
          
          if exist('fit_func','var')
           beta = nlinfit(x,y,fit_func,beta0);
           hold on;
           x2 = sort(x);
           y2=fit_func(beta,x2);
           disp('Fitting result');
           for i=1:length(beta)
           disp(beta(i));
           end
           plot(x2,real(y2),'r-');
           if any(~isreal(y2))
            plot(x2,imag(y2),'b-');
          end
           
           
          end
          
        end
    end
end


