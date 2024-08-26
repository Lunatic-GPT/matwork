function res=sliceGRAPPAKernel(S,Kx,Ky,kernel,ch)
% res=sliceGRAPPAKernel(S,Kx,Ky,kernel,ch)
% for calibration
% S: nro*npe*nsl*nch 
% the result is a nch*(length(Kx)*length(Ky)*nch)*nsl kernel
% or
% for recon
% res=sliceGRAPPAKernel(S,Kx,Ky,kernel)
% S: nro*npe*nch
% kernel: nch*(length(Kx)*length(Ky)*nch)*nsl
% res: nro*npe*nsl*nch
%Kx and Ky is a vector that specifies which neighboring k-space lines will be
%used, i.e. Kx=[-2,-1,0,1,2];
%% 2/20/2017: switch the third and forth dimensions.
S=permute(S,[1,2,4,3]);

nch=size(S,3);
if ~exist('ch','var')
    ch=1:nch;
end

nro=size(S,1);
npe=size(S,2);
Bx=length(Kx);
By=length(Ky);

if ~exist('kernel','var') || isempty(kernel)

    nsl=size(S,4);
    
S_collapse=sum(S,4);

res=zeros(length(ch),nch*Bx*By,nsl);
for i=1:nsl
    
    for j=1:length(ch)%1:nch
        
        xd= zeros(nro*npe,nch*Bx*By,'single');
        y=zeros(nro*npe,1);
        nrow=0;
        for iro=1:nro
            for ipe=1:npe
                
                
                ind1=iro+Kx;
                ind2=ipe+Ky;
                
                if ind1(1)<1  || ind2(1)<1 || ind1(end)>nro || ind2(end)>npe
                    continue;
                end
                
                
               
                nrow=nrow+1;
                xd(nrow,:)=vec(S_collapse(ind1,ind2,:));
                y(nrow)=S(iro,ipe,ch(j),i);
                
            end
        end
    
    tmp=xd(1:nrow,:)\y(1:nrow);    
    res(j,:,i)=tmp;    
    disp([i,j]);
    end
    
end

else

    nsl=size(kernel,3);
    
    res=zeros(nro,npe,length(ch),nsl);
    
    lxd=nch*By*Bx;
    
    for i=1:nsl
        tic;
        pkernel=permute(kernel(:,:,i),[2,1]);
               
        for iro=1:nro
            xd=zeros(npe,lxd);
            
            for ipe=1:npe
                
                ind1=iro+Kx;
                ind2=ipe+Ky;
                
                if ind1(1)<1  || ind2(1)<1 || ind1(end)>nro || ind2(end)>npe
                    continue;
                end
                
                
                xd(ipe,:)=reshape(S(ind1,ind2,:),[1,lxd]);
                
            end
            res(iro,:,:,i)=xd*pkernel;
        end
           
    disp(toc);          
    end
    
    res=permute(res,[1,2,4,3]);
    
end


