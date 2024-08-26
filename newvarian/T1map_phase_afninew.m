function T1map_phase_afninew(fid_prefix,ind_ex)
% T1map_sdt(fid_prefix,ind_ex[,stretch_exp])
% ind_ex: subbriks to exclude.  1 based.
% stretch_exp: default false.


if ~exist('ind_ex','var')
    ind_ex=[];
end
%epiShaper32(fid_prefix,false);

[a,info]=BrikLoad([fid_prefix,'_mag+orig']);
b=BrikLoad([fid_prefix,'_ph+orig']);

d=a.*exp(1i*b);

te2=readPar(fid_prefix,'ti');
img=readPar(fid_prefix,'image');
te2=te2(img==1);
sz= size(a);


ind_inc=setdiff(1:sz(4),ind_ex);
d_rshp = d(:,:,:,ind_inc);
te2 = te2(ind_inc);

[te2_sort,i_sort]=sort(te2);

d_rshp=d_rshp(:,:,:,i_sort);
d_rshp = 1-d_rshp(:,:,:,:)./repmat(d_rshp(:,:,:,end),[1,1,1,size(d_rshp,4)]);


t2 = zeros([sz(1:3),4]);

options=optimset('MaxIter',30,'Display','off','FunValCheck','off','TolFun',0.001,'TolX',0.001);
for i=1:size(a,1)
    disp(i);
    for j=1:size(a,2)
        for k=1:size(a,3)
            
            y = squeeze(d_rshp(i,j,k,:));
            y=abs(y);
              % [beta,r]=nlinfit(te2(:),y(:),@T1_rcvr_abs,[max(y),2,2]);
              [beta,r]=lsqcurvefit(@exp_decay2,[2,2,0],te2_sort(:),y(:),[],[],options);
               beta(4)=1;
                
            if ~any(isnan(beta)) 
             ss = sum((y-mean(y)).^2);   
             t2(i,j,k,1) = 1/beta(2);
             t2(i,j,k,4) = 1-sum(r.^2)/ss;
             t2(i,j,k,2) = beta(1);
             t2(i,j,k,3)=beta(3);
            end
            
        end
    end
end

   name = ['R1map_',fid_prefix];


write_afni(t2(:,:,:,1),name,info);

