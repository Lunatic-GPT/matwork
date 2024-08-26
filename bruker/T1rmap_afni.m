function T1rmap_afni(fid_prefix,stretch_exp)
% T2map(fid_prefix[,ind_ex,stretch_exp])
% ind_ex: subbriks to exclude.  1 based. default: [].
% stretch_exp: default false.


if ~exist('stretch_exp','var')
    stretch_exp = false;
end

if exist([fid_prefix,'_recon+orig.HEAD'],'file')   
 [a,info]=BrikLoad([fid_prefix,'_recon+orig']);
else
 [a,info]=BrikLoad([fid_prefix,'+orig']);   
end
    
te2=readbPar([fid_prefix,'/method'],'tSL');
nr=readbPar([fid_prefix,'/method'],'PVM_NRepetitions');

epi=readbPar([fid_prefix,'/method'],'Method',false);
sz= size(a);

if ~isempty(strfind(epi,'EPI'))
    nr=nr-1;
end


a=reshape(a,[sz(1:3),length(te2),nr/length(te2)]);
t2 = zeros([sz(1:3),size(a,5)]);
res=t2;
beta=t2;

options=optimset('MaxIter',30,'Display','off','FunValCheck','off','TolFun',0.001,'TolX',0.001);
for ir=1:size(a,5)
for i=1:size(a,1)
    disp(i);
    for j=1:size(a,2)
        for k=1:size(a,3)
 
            y = squeeze(a(i,j,k,:,ir));
            if ~stretch_exp
               [b,r]=lsqcurvefit(@exp_decay,[max(y),0.04],te2(:),y(:),[],[],options);
               b(3)=1;
            else
               [b,r]=lsqcurvefit(@stretch_exp_decay,[max(y),0.026,1],te2(:),y(:),[],[],options);
            end    
            if ~isnan(beta(1)) && ~isnan(beta(2)) && ~isnan(beta(3))
             ss = sum((y-mean(y)).^2);   
             t2(i,j,k,ir) = 1/b(2);
             res(i,j,k,ir) = 1-sum(r.^2)/ss;
             beta(i,j,k,ir)=b(3);
            end
            
        end
    end
end
end

write_afni(t2,[fid_prefix,'_R1r'],info);

%WriteBrikEZ(t2,info,'T1rmap',[fid_prefix,'_R1r'],'R1r~');
%WriteBrikEZ(res,info,'res_afni',[fid_prefix,'_res'],'residual~');

if stretch_exp 
  WriteBrikEZ(beta,info,'beta_afni',[fid_prefix,'_beta'],'residual~');
end

%writesdt4(t2,name);
   



