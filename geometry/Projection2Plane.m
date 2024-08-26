function vperp=Projection2Plane(v,norm)

% project vector v into a plane perpendicular to vector norm.
% v is [n1,n2,n3,..,3];
% norm is [1,3];

nd=ndims(v);
sz1=ones(1,nd);
sz1(end)=3;

norm=reshape(norm,sz1);
sz2=size(v);
sz2(end)=1;


if datenum(version('-date'))>=736580
    
vpar=norm.*sum(v.*norm,nd)./(sum(norm.^2,nd));
else
norm=repmat(norm,sz2);
vpar=norm.*repmat(sum(v.*norm,nd)./(sum(norm.^2,nd)),sz1);
end


vperp=v-vpar;



