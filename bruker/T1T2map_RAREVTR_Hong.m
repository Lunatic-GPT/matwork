function T1T2map_RAREVTR_Hong(fid_prefix,t1t2,fit_mean)
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

a=brecon2afni(fid_prefix);


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

options=optimset('MaxIter',20,'Display','off');

[tmp,ite]=min(te2);
[tmp,itr]=max(vd);
tmp=a(:,:,:,ite,itr,1);
m=tmp>0.02*max(tmp(:));

for ir=1:size(a,6)
for i=1:size(a,1)
    fprintf('%d/%d\n',i,size(a,1));
    for j=1:size(a,2)
        for k=1:size(a,3)
          if m(i,j,k)==0
              continue;
          end
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
  %     tic;
            [beta,r]=lsqcurvefit(@exp_decay2,[-max(y),t1t2_0(1),max(y)],vd(:),y(:),[],[],options);
   %        disp(toc); 
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

R1=t1;
R2=t2;


if length(vd)>=3
    
  save([fid_prefix,'_R1'],'R1');
end

if length(te2)>1
   
  save([fid_prefix,'_R2'],'R2');
end

function im=brecon2afni(d,irecon)

if ~exist('irecon','var')
    irecon=1;
end

im=read_2dseq(d,irecon);
nr=readbPar(fullfile(d,'acqp'),'NR',true);
ni=readbPar(fullfile(d,'acqp'),'NI',true);
ns=readbPar(fullfile(d,'acqp'),'NSLICES',true);
im=reshape(im,[size(im,1),size(im,2),size(im,3)/nr/ni*ns,nr*ni/ns]);


function y=exp_decay(b,x)
% y=b1*exp(-x/b2);
% %4.3f*exp(-x/%4.3f)
y=b(1)*exp(-x/b(2));

function res=readbPar(fname,par,isnum)
%res=readbPar(fname,par,isnum)

if ~exist('isnum','var')
    isnum=true;
end

fid=fopen(fname,'r');
par=['##$',par,'='];
while 1
  b=fgetl(fid);
  if b==-1
      fclose(fid);
      error([par(4:end-1), ' not found']);
  end
  
  ind=strfind(b,par);
  if ~isempty(ind) 
    res=b(ind+length(par):end);
    break;
  end

end

if res(1)=='('
    sz=str2num(res(2:end-1));
     if numel(sz)==1
            sz=[1,sz];
     end
        
    if isnum 
     res=[];
     while 1
        
      b=fgetl(fid);
      if isnum
        tmp=str2num(b);
      else
       tmp=strread(b,'%s');  
      end
      res=[res,tmp];
      
      if length(res)==prod(sz) || ~isnum
          break;
      end
      
     end
    res=reshape(res,sz(end:-1:1));
    
    else
         b=fgetl(fid);
       res=strread(b,'%s');  
     
    end
else
    if isnum    
      res=str2double(res);
    end
    
    
end
fclose(fid);
