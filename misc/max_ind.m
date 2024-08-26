function res=max_ind(data)

[tmp,ind]=max(data(:));

res=ind2subb(size(data),ind);

if length(res)==2
  res(res==1)=[];
end