function [res,pos]=distance2Line(lin,point)

% shortest distance from a point to a line defined by lin

% lin: 2*2 or 2*3 matrix that defines the two points on the line
% point [n1,n2,...,2] or [n1,n2,...,3]
% pos is the position on the line as reference to the first point.  Axis
% direction from point 1 to point 2.

sz=size(point);
point=reshape(point,[prod(sz(1:end-1)),sz(end)]);
np=size(point,1);

if datenum(version('-date'))>=736580  % later or equal to September 7, 2016; R2016b
c=point-lin(1,:);  % increased speed
else  
   c=point-repmat(lin(1,:),[np,1]);
end
b=lin(2,:)-lin(1,:);
ac=sos(c,2);
ab=sos(b,2);
% 
% tic; %%2.0
% if datenum(version('-date'))>=736580
%      an=acos(sum(c.*b,2)./ac/ab);
% else
%     an=acos(sum(c.*repmat(b,[np,1]),2)./ac/ab);
% end
% 
% res=sin(an).*ac;
% toc;


%%

%res=sqrt(1-sum(c.*b,2).^2./sum(c.^2,2)/sum(b.^2,2)).*ac;

if datenum(version('-date'))>=736580
res=sqrt(1-(sum(c.*b,2)./ac/ab).^2).*ac;
else
    res=sqrt(1-(sum(c.*repmat(b,[size(c,1),1]),2)./ac/ab).^2).*ac;
end

%%

res(isnan(res))=0;
% 
% tmp=sin(an).*ac;
% res(ac>0)=tmp(ac>0);

sz2=sz(1:end-1);
if numel(sz2)==1
    sz2=[sz2,1];
end
res=reshape(res,sz2);

if nargout>1
    if datenum(version('-date'))>=736580
      pos=sum(c.*b,2)./ac/ab.*ac;
    else
        pos=sum(c.*repmat(b,[size(c,1),1]),2)./ac/ab.*ac;
    end
pos(ac==0)=0;

pos=reshape(pos,sz2);

end










