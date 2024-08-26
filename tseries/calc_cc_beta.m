   function [cc,b,p]=calc_cc_beta(ts,ref,dim)   
   %[cc,b]=calc_cc_beta(ts,ref,dim)   
   % 
   sz=size(ts);
   
   if length(sz)==2
       sz(dim)=1;
   else
   sz(dim)=[];
   end
   if size(ts,dim)~=length(ref)
       error('dimension mismatch');
   end
   
    
    ind=1:length(sz)+1;
    ind(dim)=[];
    ind=[ind,dim];
    
    ts=permute(ts,ind);
    ts=reshape(ts,[prod(sz),length(ref)]);
    
    ref2=repmat(ref(:)'-mean(ref),[prod(sz),1]);
    

 ts=ts-repmat(mean(ts,2),[1,length(ref)]);
 
 cc=mean(ts.*ref2,2)./sqrt(mean(ts.^2,2))./sqrt(mean(ref2.^2,2));
cc=reshape(cc,sz);
if nargout>2
p=cc2tp(cc,length(ref));
end
if nargout>=2
    ref=ref(:)-mean(ref);
    ref=ref/(max(ref)-min(ref));
    x=[ref(:),ones(length(ref),1)];
    b=x\permute(ts,[2,1]);
    b=b(1,:);
    b=reshape(b,sz);
end
    