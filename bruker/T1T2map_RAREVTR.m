function T1T2map_RAREVTR(fid_prefix,t1t2,fit_mean)
% T1map_sdt(fid_prefix,[t1t2_0])
% fid_prefix.
% t1t2_0: initial values for T1 and T2. default: 2, 0.045
% fit_mean: fit averaged data; default true;
% only tested for single slice 
if ~exist('stretch_exp','var') || isempty(t1t2)
    t1t2_0 = [2,0.045];
end

if ~exist('fit_mean','var')
    fit_mean=true;
end

[a,info]=BrikLoad([fid_prefix,'+orig']);

vd=readbPar(fullfile(fid_prefix,'method'),'MultiRepetitionTime');
vd=vd/1000;
te2=readbPar(fullfile(fid_prefix,'acqp'),'ACQ_echo_time');
%vd=readbPar(fullfile(fid_prefix,'method'),'MultiRepetitionTime');
%te2=readbPar(fullfile(fid_prefix,'method'),'EffectiveTE');
NR=readbPar(fullfile(fid_prefix,'method'),'PVM_NRepetitions');

ns=readbPar(fullfile(fid_prefix,'acqp'),'NSLICES');
te2=te2/1000;
sz= size(a);

a=reshape(a,[size(a,1),size(a,2),ns,length(te2),length(vd),NR]);

if fit_mean
    a=mean(a,6);
    NR=1;
end

t2 = zeros([sz(1:2),ns,NR]);
res2=zeros([sz(1:2),ns,NR]);
t1=zeros([sz(1:2),ns,NR]);
res1=zeros([sz(1:2),ns,NR]);

yt1=mean(a,4);
yt2=mean(a,5);

options=optimset('MaxIter',30,'Display','off');

for ir=1:size(a,6)
for i=1:size(a,1)
    fprintf('%d/%d\n',i,size(a,1));
    for j=1:size(a,2)
        for k=1:size(a,3)
      
            y = squeeze(yt2(i,j,k,:,1,ir));
          
            if length(te2)>1
            [beta,r]=lsqcurvefit(@exp_decay,[max(y),t1t2_0(2)],te2(:),y(:),[],[],options);
            
            if ~any(isnan(beta)) 
             ss = sum((y-mean(y)).^2);   
             t2(i,j,k,ir) = 1/beta(2);
             res2(i,j,k,ir) = 1-sum(r.^2)/ss;
            end
            end
            
            
            y = squeeze(yt1(i,j,k,1,:,ir));
          if length(vd)>=3
              
       %       if i==54 && j==52 && k==2
            [beta,r]=lsqcurvefit(@exp_decay2,[-max(y),t1t2_0(1),max(y)],vd(:),y(:),[],[],options);
            
            if ~any(isnan(beta)) 
             ss = sum((y-mean(y)).^2);   
             t1(i,j,k,ir) = 1/beta(2);
             res1(i,j,k,ir) = 1-sum(r.^2)/ss;
            end
        %      end
          end
        end
           
    end
end

end


if length(vd)>=3
    if ~fit_mean
     write_afni(t1,[fid_prefix,'_R1'],info);
     write_afni(res1,[fid_prefix,'_R1res'],info);
    else
     write_afni(t1,[fid_prefix,'_R1_mean'],info);
     write_afni(res1,[fid_prefix,'_R1res_mean'],info);
    end 
end

if length(te2)>1
    if ~fit_mean
      write_afni(t2,[fid_prefix,'_R2'],info);
      write_afni(res2,[fid_prefix,'_R2res'],info);
    else   
      write_afni(t2,[fid_prefix,'_R2_mean'],info);
      write_afni(res2,[fid_prefix,'_R2res_mean'],info);
    end
end

