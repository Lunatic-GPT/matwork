function f=p2f(p,d1,d2)
% f=p2f(p,d1,d2)
%s = sprintf('FTest(%d,%d,x)-%f',d1,d2,p);


f=fzero(@(x) FTest(d1,d2,x)-p,13);

    