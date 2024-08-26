function [out,cm_out]=combine_over_under(under,over,cm_u,cm_o,mask)
% [out,cm_out]=combine_over_under(under,over,cm_u,cm_o,mask)
% both under and over are indices to cm_u and cm_o.
% under: n*m underlay image;  values are indices to colormap in cm_u
% over: n*m overlay image;  values are indices to colormap in cm_o
% cm_u: colormap for underlay; dimensiion ncolor*3
% cm_o: colormap for overlay; dimension ncolor*3
% mask: n*m mask for overlay; overlay color only shown for pixels whose mask value is true. 
cm_out=[cm_u;cm_o];
under(under>size(cm_u,1))=size(cm_u,1);
over(over<1)=1;

out=under;

n=size(under,1)/size(over,1);
m=size(under,2)/size(over,2);

mask=reshape(mask,[1,size(mask,1),1,size(mask,2)]);
mask=repmat(mask,[n,1,m,1]);
mask=reshape(mask,[size(mask,2)*n,size(mask,4)*m]);


over=reshape(over,[1,size(over,1),1,size(over,2)]);
over=repmat(over,[n,1,m,1]);
over=reshape(over,[size(over,2)*n,size(over,4)*m]);


out(mask>0)=over(mask>0)+size(cm_u,1);








