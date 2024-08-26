function [T2,M0]=T2map_data(data,te,mask)
% T2map_data(data,te)

if length(te)~=size(data,4)
    error('Dimension mismatch');
end

[tmp,ind]=sort(te);

im=data(:,:,:,ind);

sz=size(im);
if ~exist('mask','var')
  
   mask=ones(sz(1:3)); 
end
T2 = zeros([sz(1:3),1]);
M0=0*T2;
options=statset('FunValCheck','off');
for i=1:size(im,1)
    for j=1:size(im,2)
        for k=1:size(im,3)
            
            y = double(squeeze(im(i,j,k,:)));
          
            if mask(i,j,k)==0
                continue;
            end
            
%             if i==134 && j==161
%                 disp('');
%             else 
%                 continue;
%             end
            
         if length(te)>2
            [beta,r]=nlinfit(te(:),double(y(:)),@exp_decay,[max(y),mean(te)],options); 
            if ~isnan(beta(1)) && ~isnan(beta(2))
             %ss = sum((y-mean(y)).^2);   
             T2(i,j,k,1) = beta(2);
             M0(i,j,k,1) = beta(1);
             
 %            t2(i,j,k,2) = 1-sum(r.^2)/ss;
            end
         else
            T2(i,j,k)=1/log(y(1)/y(2))*(te(2)-te(1));
            
            M0(i,j,k) = y(1)*exp(te(1)/T2(i,j,k));
            
         end
        end
    end
end

%write_afni(t2(:,:,1,1),name);
   



