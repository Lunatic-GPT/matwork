function omat=combine_flirt_omat(f1,f2)
% combining the transformation matrices of f1 and f2, f1 was applied first

if isa(f1,'char')
omat1=load(f1);
omat1=omat1(1:3,:);
end


if isa(f2,'char')
omat2=load(f2);
omat2=omat2(1:3,:);
end


omat=omat2*omat1;