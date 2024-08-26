function res = adjD(y)


res = adjDx(y(:,:,:,:,1)) + adjDy(y(:,:,:,:,2));

if size(y,3)>1
  %  res=res+adjDz(y(:,:,:,:,3));
end



function res = adjDy(x)
res = x(:,[end,1:end-1],:,:) - x;

function res = adjDx(x)
res = x([end,1:end-1],:,:,:) - x;


function res = adjDz(x)
res = x(:,:,[end,1:end-1],:) - x;
