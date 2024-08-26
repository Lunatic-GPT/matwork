function [T1,resid]=T1map(a,TI)
% T1map(fid_prefix,ind_ex[,stretch_exp])
% ind_ex: subbriks to exclude.  1 based. default:[].
% stretch_exp: default false.

sz=size(a);

T1 = zeros(sz(1:3));
resid=T1;

options=statset('FunValCheck','off');

for i=1:size(a,1)
    for j=1:size(a,2)
        for k=1:size(a,3)
            
            y = squeeze(a(i,j,k,:));
          
            if ~any(y>0)
                continue;
            end
            
               [beta,r]=nlinfit(TI(:),y(:),@exp_decay2,[-max(y),mean(TI),max(y)],options);
               
            if ~any(isnan(beta)) 
             resid(i,j,k) = rms(exp_decay2(beta,TI(:))-y(:));   
             
             T1(i,j,k)=beta(2);
             
            end
            
        end
    end
    disp(i);
end





