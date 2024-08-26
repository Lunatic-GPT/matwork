function cstable_epi_mTR(N,nTRcs,nc,reduction,ns,suffix)
% cstable_epi_mTR(N,nTRcs,nc,reduction,ns,suffix)
%N=128;
%nTRcs=300;
%nc=6;  %includes the edge line
%reduction = 2;
%ns=1;

nout =(N/reduction-nc);

blips=[];

oint = [2:N/2-nc/2,N/2+nc/2:N];


blips_all = zeros((N/reduction-1),ns,nTRcs);


mask_all = zeros((N/reduction-1),ns,nTRcs);

for is=1:ns

   j=1; 
  while j<=nTRcs/reduction
   k=randperm(length(oint));
  
    for i=1:reduction
      
        ind=[1,N/2-nc/2+1:N/2+nc/2-1];
        tmp=k((i-1)*nout+1:i*nout);
        ind(end+1:end+nout)=oint(tmp);
        ind=sort(ind);
        
        mask = zeros(N,N);
        mask(:,ind)=1;
   
        if any(diff(ind)>10)       
            break;
        end
       
        blips_all(:,is,(j-1)*reduction+i)=diff(ind);  
         mask_all(ind,is,(j-1)*reduction+i)=1;  
        
        disp([is,j,i]);
    end
    
     if any(diff(ind)>10)       
         continue;
     else       
       j=j+1;
     end
  end
  figure;imshow(squeeze(mask_all(:,is,:)),[]);
  %
end

%%
fname = sprintf('blipfactor_R%d_N%d_nTR%d_ns%d_nc%d_%s',reduction,N,nTRcs,ns,nc,suffix);


if exist(fname,'file')
    warning('file already exist');
else
fid=fopen(fname,'w');


blips_all=blips_all(:);
for i=1:length(blips_all)
  fprintf(fid,'%d\n',blips_all(i));
end
fclose(fid);



end