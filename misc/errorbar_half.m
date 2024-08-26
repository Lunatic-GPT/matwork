function errorbar_half(x,y,e,w,tb,varargin)
%errorbar_half(x,y,e,w,tb,varargin)
hold on;

for i=1:length(x(:))
    
    sym=varargin{1};
    if length(sym)==3 || sym(2)~='-'
        sym=[sym(1),'-'];
    end
    varargin{1}=sym;
    
 if strcmp(tb,'top') || strcmp(tb,'both')   
     k=1;
    plot([x(i),x(i)],[y(i),y(i)+e(i)*k],varargin{:});
    plot([x(i)-w/2,x(i)+w/2],[y(i)+e(i)*k,y(i)+e(i)*k],varargin{:});
 end

  if strcmp(tb,'bottom') || strcmp(tb,'both')   
     k=-1;
    plot([x(i),x(i)],[y(i),y(i)+e(i)*k],varargin{:});
    plot([x(i)-w/2,x(i)+w/2],[y(i)+e(i)*k,y(i)+e(i)*k],varargin{:});
 end

end

