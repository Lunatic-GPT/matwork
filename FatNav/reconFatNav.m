function im=reconFatNav(data,par)
% center give the center position in the readout (I2S), PE (A2P), and
% Partition (R2L)
% directions in units of the fov in those directions.


Lin=par.Lin; %0 based
Par=par.Par; %0 based
Rep=par.Rep; %0 based
iRep=par.iRep; %1 based
iRep_ref=par.iRep_ref; %1 based
center=par.center; %center along readout, pe, and par;
prefix=par.prefix;
maxLin=par.maxLin; %0 based; without undersampling and partial Fourier 
maxPar=par.maxPar; %0 based

mref=zeros(maxLin+1,maxPar+1);
nro=size(data,1);
ndim=ndims(data);


if ndim==2
    Lin=Lin(1:32:end);
    Par=Par(1:32:end);
    Rep=Rep(1:32:end);
    data=reshape(data,[size(data,1),32,length(Lin)]);
end


d2=zeros(nro,maxLin+1,maxPar+1,length(iRep),32,'single');
dref=zeros(nro,maxLin+1,maxPar+1,32);

for i=1:length(Lin)

    ind=find(Rep(i)+1==iRep);
    if ~isempty(ind)  
        % the sign inside exp is determined as 
        % 1.  shift the object in the negative direction -center direction
        % 2. the axes for the recon image and the axes for center are
        % opposite; another - sign
        % 3. another - sign to shift data x, i.e. -ikx
                                       
     d2(:,Lin(i)+1,Par(i)+1,ind,:)= data(:,:,i)*exp(-1i*2*pi*(center(3)*double(Par(i))+center(2)*double(Lin(i)))); 
        
    end
    
    if (Rep(i)+1==iRep_ref)
       dref(:,Lin(i)+1,Par(i)+1,:)= data(:,:,i)*exp(-1i*2*pi*(center(3)*double(Par(i))+center(2)*double(Lin(i))));    
       mref(Lin(i)+1,Par(i)+1)=mref(Lin(i)+1,Par(i)+1)+1;
    end
    
   
end


m2=clusterize2(mref,10);
s2=sum(m2,2);
s1=sum(m2,1);

acs_row = s2>=max(s2)-2;
acs_col=s1>=max(s1)-2;

irow=1;
while 1
    if any(mref(irow,:)>0)
        break;
    end
    irow=irow+1;
end
inc_row=irow:size(mref,1);

icol=1;
while 1
    if any(mref(:,icol)>0)
        break;
    end
    icol=icol+1;
end
inc_col=icol:size(mref,2);

dacs=dref(:,acs_row,acs_col,:);


im=zeros(nro,maxLin+1,maxPar+1,length(iRep),'single');

%iRep=;     
tic;
    for j=1:nro
        tmp=shiftdim(d2(j,:,:,:,:),1);    
 
        tmp=permute(tmp,[1,2,4,3]);
        kres=0*tmp;
        %kres=GRAPPA(permute(tmp,[1,2,4,3]),squeeze(dacs(j,:,:,1,:)),[5,5],0.01,[4,4]);
        fprintf('%d/%d:',j,nro);
        
  
       %tmp2=GRAPPA(tmp(inc_row,inc_col,:,:),squeeze(dacs(j,:,:,:)),[7,7],0.01);
      
        %   tmp1=GRAPPA_new(tmp(inc_row,inc_col,:,:),squeeze(dacs(j,:,:,:)),[7,7],0.01,[4,4]);
         
             kres(inc_row,inc_col,:,:)=GRAPPA_new(tmp(inc_row,inc_col,:,:),squeeze(dacs(j,:,:,:)),[7,7],0.01,[4,4]);
       %   kres(inc_row,inc_col,:,:)=GRAPPA(tmp(inc_row,inc_col,:,:),squeeze(dacs(j,:,:,:)),[7,7],0.01);
      
       
        im(j,:,:,:) = sos(ifft1c(ifft1c(kres,1),2),3); %need to check before
        time_left(j,nro,toc);
    end

if nargout==0
    d=im;  %rename to d
    save([prefix,'.mat'],'d','iRep');
end




    
    