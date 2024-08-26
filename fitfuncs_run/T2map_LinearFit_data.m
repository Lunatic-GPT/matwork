function [T2,M0]=T2map_LinearFit_data(data,te,mask)
% T2map_LinearFit_data(data,te[,mask])
% ind_ex: subbriks to exclude.  1 based. default: [].
% stretch_exp: default false.

if length(te)~=size(data,4)
    error('Dimension mismatch');
end

[te,ind]=sort(te);

im=data(:,:,:,ind);

sz=size(im);
if ~exist('mask','var')
  
   mask=ones(sz(1:3)); 
end
T2 = zeros([sz(1:3),1]);
M0 = zeros([sz(1:3),1]);

for i=1:size(im,1)
     disp(i);
    for j=1:size(im,2)
       
        for k=1:size(im,3)
            
            
            if ~(i==244 && j==250  && k==14)
                
            %    continue;
            end
            
            if mask(i,j,k)==0
                continue;
            end
            y = double(squeeze(im(i,j,k,:)));
            
            if any(y==0)
                continue;
            end
            
            logy=log(y);
            
            x=cat(2,ones(length(te),1),-te(:));
            
            y2=logy.*y;
            
            x2=x.*repmat(y,[1,2]);
            
            b=x2\y2;
            
            T2(i,j,k)=1/b(2);
            M0(i,j,k)=exp(b(1));
            
            %% calculate the weights for the data
            
            
            
        end
        
        
    end
end

   



