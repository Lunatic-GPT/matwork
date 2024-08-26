function abc=tlrc2orig(x,y,z,f_tlrc)
%abc=tlrc2orig(x,y,z,f_tlrc)

cmd = sprintf('cat_matvec %s::WARP_DATA',f_tlrc);
[s,m]=unix(cmd);
if s==0
m = str2num(m);
abc = inv(m(:,1:3))*([x;y;z]-m(:,4));
%abc=m(:,1:3)*[x;y;z]+m(:,4);
else
    disp(m);
end