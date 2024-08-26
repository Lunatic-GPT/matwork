function b1=B1map_th2th(d,theta)
%b1=B1map_th2th(d,theta)
% dlist should contain two file names which corresponds to theta and 2 theta, respectively.
% or a matrix with the fourth dimension = 2




dmax=max(abs(d(:)));
    b1=zeros(size(d(:,:,:,1)));
    
for i=1:size(d,1)
    for j=1:size(d,2)
        for k=1:size(d,3)
            y=d(i,j,k,:);
            y=y(:);
            
           
            if i==107 && j==61
                disp('');
            end
           %if max(abs(y))<dmax*0.03
           if abs(y(1))==0
               disp([i,j]);
                 continue;
           end
            
                
            
            r=abs(y(2)/y(1));
            if abs(angle(y(2)/y(1)))>pi/2    
                r=-abs(y(2)/y(1));
            end
            x=(r+sqrt(r^2+8))/4;
            if abs(x)>1 
                continue;
            end
            b1(i,j,k)=acos(x)*90/theta*180/pi;
            
        end
    end
end

b1=abs(b1);
save(sprintf('B1map_th2th_%d',theta),'b1');

