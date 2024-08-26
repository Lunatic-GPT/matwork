function [R2,b,Rsquare]=T2_ROI(fid_prefix,m,ind_ex,stretch_exp)
% [t2,b,Rsquare]=T2_ROI(fid_prefix,m,ind_ex,stretch_exp)
% ind_ex: subbriks to exclude.  1 based. default: [].
% stretch_exp: default false.

if ~exist('stretch_exp','var')
    stretch_exp = false;
end

if ~exist('ind_ex','var')
    ind_ex = [];
end

a=mean_roi([fid_prefix,'+orig'],m);
ref_pos=readPar([fid_prefix,'.fid'],'ref_pos');

if ~any(ind_ex==ref_pos+1)
    warning('The reference subbrik %d is also excluded',ref_pos+1);
    ind_ex=[ind_ex,ref_pos+1];
end


te2=parValArray(fid_prefix,'te');

arraydim = readPar(fid_prefix,'arraydim');
if length(a)==arraydim-1
   b = cat(4,a(:,:,:,1:ref_pos+1),a(:,:,:,ref_pos+1:end));
elseif length(a)~=arraydim
  error('Data set error');
else
    b=a;
end

ind_inc=setdiff(1:length(a),ind_ex);
y = b(ind_inc);
te2 = te2(ind_inc);

R2=0;
Rsquare=0;
b=0;
options=statset('FunValCheck','off');
            
            if ~stretch_exp
               [beta,r]=nlinfit(te2(:),y(:),@exp_decay,[max(y),0.026],options);
               beta(3)=1;
            else
               [beta,r]=nlinfit(te2(:),y(:),@stretch_exp_decay,[max(y),0.026,1],options);
            end 
            
            if ~isnan(beta(1)) && ~isnan(beta(2)) && ~isnan(beta(3))
             ss = sum((y-mean(y)).^2);   
             R2 = 1/beta(2);
             Rsquare = 1-sum(r.^2)/ss;
             b=beta(3);
             
        %    figure;
         %   semilogy(te2,y,'o');
            
         %   hold on;
         %   x=linspace(0,max(te2));
         %   semilogy(x,stretch_exp_decay(beta,x),'r-');
            end



