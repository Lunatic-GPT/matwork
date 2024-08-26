function T2map_afni(fid_prefix,stretch_exp)
% T2map(fid_prefix[,ind_ex,stretch_exp])
% ind_ex: subbriks to exclude.  1 based. default: [].
% stretch_exp: default false.

if ~exist('stretch_exp','var')
    stretch_exp = false;
end

[a,info]=BrikLoad([fid_prefix,'+orig']);

te2=readPar(fid_prefix,'te');
img=readPar(fid_prefix,'image');
te2=te2(img==1);    
sz= size(a);
t2 = zeros([sz(1:3),3]);

options=optimset('MaxIter',30,'Display','off','FunValCheck','off','TolFun',0.001,'TolX',0.001);
for i=1:size(a,1)
    disp(i);
    for j=1:size(a,2)
        for k=1:size(a,3)
            
            y = squeeze(a(i,j,k,:));
          
            if ~stretch_exp
               [beta,r]=lsqcurvefit(@exp_decay,[max(y),0.026],te2(:),y(:),[],[],options);
               beta(3)=1;
            else
               [beta,r]=lsqcurvefit(@stretch_exp_decay,[max(y),0.026,1],te2(:),y(:),[],[],options);
            end    
            if ~isnan(beta(1)) && ~isnan(beta(2)) && ~isnan(beta(3))
             ss = sum((y-mean(y)).^2);   
             t2(i,j,k,1) = 1/beta(2);
             t2(i,j,k,2) = 1-sum(r.^2)/ss;
             t2(i,j,k,3)=beta(3);
            end
            
        end
    end
end


if ~stretch_exp
   name = ['R2map_',fid_prefix];
else    
   name = ['R2mapStretch_',fid_prefix];
end

WriteBrikEZ(t2,info,'T2map_afni',name,'R2~R2~beta~');
%writesdt4(t2,name);
   



