function dsb2=reorder_fl_fq_data_interp(dsb,istart,Set,nphase,method)
% dsb should be in the format of [nro*npe*nsl*nch]; 
% peaks
% the index for starting a new k-space line
% output: [nro,lPhaseEncodingLines,nsl,nch,nvenc,nphase];

if ~exist('method','var')
    method='linear';
end

istart=[1,istart(:)'];

nro=size(dsb,1);

nk=length(istart)-1;
nvenc=length(unique(Set));
ns=size(dsb,3);
nch=size(dsb,4);

dsb2=zeros(nphase,nro,nk,ns,nch,nvenc,'single');

y=permute(dsb,[2,1,3,4]);
xi=linspace(0,1,nphase+1);
    
for i=1:length(istart)-1
    for j=1:nvenc
    disp(i);
    
    
      
      s=Set(istart(i):istart(i+1)-1);

      nt=sum(s==j-1);
      x=linspace(0,1,nt+1);
      x=x(1:end-1);
      
      ytmp=y(istart(i):istart(i+1)-1,:,:,:);
      
      dsb2(:,:,i,:,:,j)=interp1(x,ytmp(s==j-1,:,:,:),xi(1:end-1),method,0);
    
    end
end

dsb2=permute(dsb2,[2,3,4,5,6,1]);









  
  
