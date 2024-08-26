function d=setv_roi(d,roi,val)
% ndims(d)>ndims(roi); the first dimensioins should match.
% val should match the size of d for the remaining dimensions.
nd=ndims(roi);
nvox=sum(roi(:)>0);
sz=size(d);
sz(1:nd)=1;
roi=repmat(roi,sz);

val=shiftdim(val,-1);

val=repmat(val,[nvox,ones(1,ndims(val))]);
d(roi>0)=val;



