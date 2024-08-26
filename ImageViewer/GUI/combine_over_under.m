function [out,cm_out]=combine_over_under(under,over,cm_u,cm_o,mask)
% [out,cm_out]=combine_over_under(under,over,cm_u,cm_o,mask)
% both under and over are indices to cm_u and cm_o.
% under: n*m underlay image;  values are indices to colormap in cm_u
% over: n*m overlay image;  values are indices to colormap in cm_o
% cm_u: colormap for underlay; dimensiion ncolor*3
% cm_o: colormap for overlay; dimension ncolor*3
% mask: n*m mask for overlay; overlay color only shown for pixels whose mask value is true. 
% for 
% Note: In imshow, the colormap index is 0 based when uint16 but 1 based when
% single or double; 
% Here out will be double.
if isempty(cm_u)
    cm_u=gray(100);
    under=under/max(under(:))*100;
end

cm_out=[cm_u;cm_o];

if isa(under,'double') || isa(under,'single')
under(under>size(cm_u,1))=size(cm_u,1);
else
under(under>size(cm_u,1)-1)=size(cm_u,1)-1;
under=under+1;
under=double(under);
end

over(over<1)=1;

n=size(under,1)/size(over,1);
m=size(under,2)/size(over,2);

out=under;
for i=1:size(out,3)
mask2=reshape(mask(:,:,i),[1,size(mask,1),1,size(mask,2)]);
mask2=repmat(mask2,[n,1,m,1]);
mask2=reshape(mask2,[size(mask2,2)*n,size(mask2,4)*m]);

over2=reshape(over(:,:,i),[1,size(over,1),1,size(over,2)]);
over2=repmat(over2,[n,1,m,1]);
over2=reshape(over2,[size(over2,2)*n,size(over2,4)*m]);

out2=out(:,:,i);
out2(mask2>0)=over2(mask2>0)+size(cm_u,1);

out(:,:,i)=out2;

end






