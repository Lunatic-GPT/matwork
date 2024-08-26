function  CBF_SVD(fname,bl,Ca,thr,outfname,fname_bl)
% CBF_SVD( prefix,bl,Ca,thr,suffix/or voxel index)
% fname: file name or data
% bl: baseline time points
% Ca: arterial input function
% thr: threshold for non-zero components relative to the first component.
% outfname: output file name or voxel index for testing
% 
% This program implments the SVD method in:
% Leif Ostergaard et al., MRM 36715-725 (1996)

tic;
if isa(Ca,'char')
    Ca=load(Ca);
end

if ~isempty(outfname) && isa(outfname,'double')
    debug=true;
    idebug=outfname;
else
    debug=false;
end

%Ca=1-Ca/max(Ca);
Ca=Ca-mean(Ca(bl));  
[tmp,ind_pk]=max(abs(Ca));
if Ca(ind_pk)<0
    Ca = -Ca;
end
Ca=Ca/abs(Ca(ind_pk));


if isa(fname,'char')
d=ri(fname);
else
    d=fname;
end

d=double(d);
sz=size(d);
if exist('fname_bl','var')
    d_bl=ri(fname_bl);
    
    bl=mean(d_bl(:,:,:,bl),4);
    bl=repmat(bl,[1,1,1,sz(4)]);
    
    d=-d./bl;
    
else
    bl=mean(d(:,:,:,bl),4);
    bl=repmat(bl,[1,1,1,sz(4)]);
    
    d=1-d./bl;
    
end

nt=size(d,4);
a=zeros(nt,nt);
if ~debug
%matlabpool(4);
end

for i=1:size(d,4)
    for j=1:size(d,4)     
      if j<=i
          s=zeros(1,3);
          n=4;
          if (i-j)>0
           s(1)=Ca(i-j);
           n=n+1;
          end
          if i-j+2<=nt
              s(3)=Ca(i-j+2);
              n=n+1;
          end
          s(2)=Ca(i-j+1);
          a(i,j)=(s(1)+s(3)+4*s(2))/n;
      end
      
    end
end

[u,s,v]=svd(a);

s=diag(s);
if thr<1
thr=abs(s(1))*thr;
s(abs(s)<thr)=0;

fprintf('non-zero components = %d\n',sum(abs(s)>=thr));
%s(1:end-nnz)=0;
is=1./s;
is(abs(s)<thr)=0;
else
    nnz=thr;
    is=1./s;
    s(nnz+1:end)=0;
    is(nnz+1:end)=0;
end

is=diag(is);
s=diag(s);

h=zeros(size(d));
fit=zeros(size(d));
for i=1:size(d,1)
    if ~debug
      disp(i);
    end
    
    tmp=zeros(size(d,2),size(d,3),size(d,4));
    
    fittmp=zeros(size(d,2),size(d,3),size(d,4));
    for j=1:size(d,2)
        for k=1:size(d,3)
            if debug && any(idebug~=[i,j,k])
               continue;
            end
          y=squeeze(d(i,j,k,:));
          tmp(j,k,:)=v*is*u'*y;
          fittmp(j,k,:)=u*s*v'*squeeze(tmp(j,k,:));                    
        end
    end
    
    h(i,:,:,:)=tmp;
    fit(i,:,:,:)=fittmp;
    
end

if ~debug
%  matlabpool close;
end

%ind_inc=setdiff(1:size(fit,4),bl);
ind_inc=1:size(fit,4);
CBV=mean(fit(:,:,:,ind_inc),4)./mean(Ca(ind_inc));
%CBV=mean(d(:,:,:,ind_inc),4)/mean(Ca(ind_inc));

CBF=max(h,[],4);
MTT=CBV./CBF;  %% in units of TR

if debug
    figure;
    
    subplot(1,3,1);
    plot(Ca);
    set(gca,'FontSize',14);
    title('AIF');
    xlabel('Time (TR)');
    
    subplot(1,3,2);plot(squeeze(d(idebug(1),idebug(2),idebug(3),:)));
    hold on;
    plot(squeeze(fit(idebug(1),idebug(2),idebug(3),:)),'r');
    set(gca,'FontSize',14);
    title('Time Course');
    
    legend('Measured','SVD Fit');
    ylabel('Intensity Change');
    xlabel('Time (TR)');
    
    subplot(1,3,3);plot(squeeze(h(idebug(1),idebug(2),idebug(3),:)));
    set(gca,'FontSize',14);
    
    title('Residual Function');
    xlabel('Time (TR)');
 %  ylabel('Intensity Change');
    set(gca,'FontSize',14);
    fprintf('CBF = %f; CBV = %f; MTT = %f\n',CBF(idebug(1),idebug(2),idebug(3)),CBV(idebug(1),idebug(2),idebug(3)),...
        MTT(idebug(1),idebug(2),idebug(3)));
    return;
end


    save([outfname,'_irf'],'h');
    save([outfname,'_fit'],'fit');
    save([outfname,'_CBV'],'CBV');
    save([outfname,'_MTT'],'MTT');
    save([outfname,'_CBF'],'CBF');


fprintf('CBF_SVD finished in %4.3f s\n',toc);
