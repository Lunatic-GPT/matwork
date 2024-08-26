function T1map_sdtnew(fid_prefix,ind_ex,stretch_exp)
% T1map_sdt(fid_prefix,ind_ex[,stretch_exp])
% ind_ex: subbriks to exclude.  1 based.
% stretch_exp: default false.

if ~exist('stretch_exp','var')
    stretch_exp = false;
end

if ~exist('ind_ex','var')
    ind_ex=[];
end

a=rdSdt(fid_prefix);

%te2=readPar(fid_prefix,'ti');
%img=readPar(fid_prefix,'image');
%te2=te2(img==1);
te2=readPar(fid_prefix,'transit');
sz= size(a);


b=a;
ind_inc=setdiff(1:sz(4),ind_ex);
b_rshp = b(:,:,:,ind_inc);
te2 = te2(ind_inc);

t2 = zeros([sz(1:3),4]);
options=optimset('MaxIter',30,'Display','off');
for i=1:size(a,1)
    disp(i);
    for j=1:size(a,2)
        for k=1:size(a,3)
            
            y = squeeze(b_rshp(i,j,k,:));
          
            if ~stretch_exp
              % [beta,r]=nlinfit(te2(:),y(:),@T1_rcvr_abs,[max(y),2,2]);
              [beta,r]=lsqcurvefit(@T1_rcvr_abs,[max(y),2,2],te2(:),y(:),[],[],options);
               beta(4)=1;
            else
              % [beta,r]=nlinfit(te2(:),y(:),@T1_rcvr_noexp,[max(y),2,2,1]);
            end    
            if ~any(isnan(beta)) 
             ss = sum((y-mean(y)).^2);   
             t2(i,j,k,1) = 1/beta(3);
             t2(i,j,k,2) = 1-sum(r.^2)/ss;
             t2(i,j,k,3) = beta(2);
             t2(i,j,k,4)=beta(4);
            end
            
        end
    end
end
if ~stretch_exp
   name = ['R1map_',fid_prefix];
else    
   name = ['R1mapStretch_',fid_prefix];
end


writesdt4(t2(:,:,:,1),name);

