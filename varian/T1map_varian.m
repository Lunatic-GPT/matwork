function T1map(fid_prefix,ind_ex,stretch_exp)
% T1map(fid_prefix,ind_ex[,stretch_exp])
% ind_ex: subbriks to exclude.  1 based. default:[].
% stretch_exp: default false.

if ~exist('stretch_exp','var')
    stretch_exp = false;
end

if ~exist('ind_ex','var')
    ind_ex=[];
end


[a,info]=BrikLoad([fid_prefix,'_mag+orig']);
ref_pos=readPar([fid_prefix,'.fid'],'ref_pos');

if ~any(ind_ex==ref_pos+1)
    warning('The reference subbrik %d is also excluded',ref_pos+1);
    ind_ex=[ind_ex,ref_pos+1];
end

te2=parValArray(fid_prefix,'ti_2');
    
sz= size(a);

arraydim = readPar(fid_prefix,'arraydim');
if size(a,4)==arraydim-1
   b = cat(4,a(:,:,:,1:ref_pos+1),a(:,:,:,ref_pos+1:end));
elseif size(a,4)~=arraydim
  error('Data set error');
else
    b=a;
end
ind_inc=setdiff(1:sz(4),ind_ex);
b_rshp = b(:,:,:,ind_inc);
te2 = te2(ind_inc);

t2 = zeros([sz(1:3),4]);
options=statset('FunValCheck','off');
for i=1:size(a,1)
    for j=1:size(a,2)
        for k=1:size(a,3)
            
            y = squeeze(b_rshp(i,j,k,:));
          
            if ~stretch_exp
               [beta,r]=nlinfit(te2(:),y(:),@T1_rcvr_abs,[max(y),2,2],options);
               beta(4)=1;
            else
               [beta,r]=nlinfit(te2(:),y(:),@T1_rcvr_noexp,[max(y),2,2,1],options);
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

labels = 'R1~Rsquare~Mint~beta~';
history = sprintf('T1map(%s,ind_ex=%s,stretch_exp=%d)',fid_prefix,num2str(ind_ex),stretch_exp);

if ~stretch_exp
   name = ['R1map_abs_',fid_prefix];
else    
   name = ['R1mapStretch_',fid_prefix];
end

WriteBrikEZ(t2,info,history,name,labels);
   



