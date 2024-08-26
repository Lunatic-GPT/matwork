function res=calc_FA(l)

lm=mean(l,4);

nom=sqrt((l(:,:,:,1)-lm).^2+(l(:,:,:,2)-lm).^2+(l(:,:,:,3)-lm).^2);
den=sqrt(l(:,:,:,1).^2+l(:,:,:,2).^2+l(:,:,:,3).^2);

res=sqrt(3/2)*nom./den;

