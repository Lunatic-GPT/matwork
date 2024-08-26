function x = fnlCg_matlab(x0,params)
%-----------------------------------------------------------------------
%
% res = fnlCg(x0,params)
%
% implementation of a L1 penalized non linear conjugate gradient reconstruction
%
% The function solves the following problem:
%
% given k-space measurments y, and a fourier operator F the function 
% finds the image x that minimizes:
%
% Phi(x) = ||F* W' *x - y||^2 + lambda1*|x|_1 + lambda2*TV(W'*x) 
%
%
% the optimization method used is non linear conjugate gradient with fast&cheap backtracking
% line-search.
% 
% (c) Michael Lustig 2007
%-------------------------------------------------------------------------
tic;
x = x0;


% line search parameters
maxlsiter = params.lineSearchItnlim ;
gradToll = params.gradToll ;
alpha = params.lineSearchAlpha;    beta = params.lineSearchBeta;
t0 = params.lineSearchT0;
k = 0;

% copmute g0  = grad(Phi(x))
XFM=params.XFM;
TV=params.TV;
FT=params.FT;
TVWeight=params.TVWeight;
xfmWeight=params.xfmWeight;
p=params.pNorm;
data=params.data;
l1Smooth=params.l1Smooth;

g0=wGradient(x,XFM,FT,TV,data,xfmWeight,TVWeight,p,l1Smooth);

dx = -g0;
% iterations



% backtracking line-search

	% pre-calculate values, such that it would be cheap to compute the objective
	% many times for efficient line-search

	[FTXFMtx, FTXFMtdx, DXFMtx, DXFMtdx] = preobjective2(x, dx, XFM,TV,FT,TVWeight);%preobjective(x, dx, params);
   % [FTXFMtx, FTXFMtdx, DXFMtx, DXFMtdx] = preobjective(x, dx, params);

    fmin  =@(x) objective2(FTXFMtx, FTXFMtdx, DXFMtx, DXFMtdx, x,dx,0, p,data,TVWeight,xfmWeight,l1Smooth);    
    
    options=optimset('MaxIter',50,'Display','iter','UseParallel','always','TolFun',0.0001,'Jacobian','off','Algorithm','active-set');
     
    x=fminunc(fmin,x0,options);
    
          


function [FTXFMtx, FTXFMtdx, DXFMtx, DXFMtdx] = preobjective(x, dx, params)

% precalculates transforms to make line search cheap

FTXFMtx = params.FT*(params.XFM'*x);
FTXFMtdx = params.FT*(params.XFM'*dx);

if params.TVWeight
    DXFMtx = params.TV*(params.XFM'*x);
    DXFMtdx = params.TV*(params.XFM'*dx);
else
    DXFMtx = 0;
    DXFMtdx = 0;
end


function [FTXFMtx, FTXFMtdx, DXFMtx, DXFMtdx] = preobjective2(x, dx, XFM,TV,FT,TVWeight)

% precalculates transforms to make line search cheap

FTXFMtx = FT*(XFM'*x);
FTXFMtdx = FT*(XFM'*dx);

if TVWeight
    DXFMtx = TV*(XFM'*x);
    DXFMtdx = TV*(XFM'*dx);
else
    DXFMtx = 0;
    DXFMtdx = 0;
end


function [res, obj, RMS,TV,XFM] = objective(FTXFMtx, FTXFMtdx, DXFMtx, DXFMtdx, x,dx,t, params)
%calculated the objective function

p = params.pNorm;

obj = FTXFMtx + t*FTXFMtdx - params.data;
obj = obj(:)'*obj(:);

if params.TVWeight
    w = DXFMtx(:) + t*DXFMtdx(:);
    TV = (w.*conj(w)+params.l1Smooth).^(p/2); 
else
    TV = 0;
end

if params.xfmWeight
   w = x(:) + t*dx(:); 
   XFM = (w.*conj(w)+params.l1Smooth).^(p/2);
else
    XFM=0;
end

TV = sum(TV.*params.TVWeight(:));
XFM = sum(XFM.*params.xfmWeight(:));
RMS = sqrt(obj/sum(abs(params.data(:))>0));

res = obj + (TV) + (XFM) ;

function [res, obj, RMS,TV,XFM] = objective2(FTXFMtx, FTXFMtdx, DXFMtx, DXFMtdx, x,dx,t, p,data,TVWeight,xfmWeight,l1Smooth)
%calculated the objective function

%p = params.pNorm;

obj = FTXFMtx + t*FTXFMtdx - data;
obj = obj(:)'*obj(:);

if TVWeight
    w = DXFMtx(:) + t*DXFMtdx(:);
    TV = (w.*conj(w)+l1Smooth).^(p/2); 
else
    TV = 0;
end

if xfmWeight
   w = x(:) + t*dx(:); 
   XFM = (w.*conj(w)+l1Smooth).^(p/2);
else
    XFM=0;
end


TV = sum(TV.*TVWeight(:));
XFM = sum(XFM.*xfmWeight(:));
RMS = sqrt(obj/sum(abs(data(:))>0));

res = obj + (TV) + (XFM) ;

function grad = wGradient(x,XFM,FT,TV,data,xfmWeight,TVWeight,p,l1Smooth)

gradXFM = 0;
gradTV = 0;

gradObj = gOBJ2(x,XFM,FT,data);
if xfmWeight
gradXFM = gXFM2(x,p,l1Smooth);
end
if TVWeight
 gradTV=gTV2(x,p,TV,XFM,l1Smooth);
end

grad = (gradObj +  xfmWeight.*gradXFM + TVWeight.*gradTV);


function gradObj = gOBJ2(x,XFM,FT,data)
% computes the gradient of the data consistency

	gradObj = XFM*(FT'*(FT*(XFM'*x) - data));

gradObj = 2*gradObj ;


function gradObj = gOBJ(x,params)
% computes the gradient of the data consistency

	gradObj = params.XFM*(params.FT'*(params.FT*(params.XFM'*x) - params.data));

gradObj = 2*gradObj ;

function grad = gXFM(x,params)
% compute gradient of the L1 transform operator

p = params.pNorm;

grad = p*x.*(x.*conj(x)+params.l1Smooth).^(p/2-1);


function grad = gXFM2(x,p,l1Smooth)
% compute gradient of the L1 transform operator

grad = p*x.*(x.*conj(x)+l1Smooth).^(p/2-1);


function grad = gTV(x,params)
% compute gradient of TV operator

p = params.pNorm;

Dx = params.TV*(params.XFM'*x);

G = p*Dx.*(Dx.*conj(Dx) + params.l1Smooth).^(p/2-1);
grad = params.XFM*(params.TV'*G);


function grad = gTV2(x,p,TV,XFM,l1Smooth)
% compute gradient of TV operator

Dx = TV*(XFM'*x);

G = p*Dx.*(Dx.*conj(Dx) + l1Smooth).^(p/2-1);
grad = XFM*(TV'*G);





