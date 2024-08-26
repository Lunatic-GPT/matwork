function res=getv_roi(d,roi)
% ndims(d)>ndims(roi); the first dimensioins should match.
% output: nvox*size_of_remaining_dim
nd=ndims(roi);
sz=size(d);
sz(1:nd)=1;
nvox=sum(roi(:)>0);
roi=repmat(roi,sz);
res=d(roi>0);

sz(1)=nvox;
res=squeeze(reshape(res,sz));


