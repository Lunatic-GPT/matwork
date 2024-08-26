function s_tissue=TotalSignal_Tissue_PerpPlane(a,dchi,B0,TE,theta,R)
% B0: unit T
% TE: unit ms
% rhoc is the density at center normalized by rho in tissue; density
% outside is assumed to be 1.
% theta: degree
% R: radius of the ROI
% a: radius of the vessels
% dchi: ppm % changed 7/19/2017
% output s: also weighted by rhoc; The signal is an average over a circle with
% diameter R and perpendicular to the vessel.
% results proportional to R^2 and is a function of a/R.
% Results from: Hsieh; MRI, 33:420-436 (2015)

gamma=42.58e6*2*pi;  %rad/T/s
theta=theta/180*pi;

TE = TE/1000;

dchi=dchi*1e-6;

g=0.5*gamma*dchi*B0*TE;

gp=g*sin(theta);
s_tissue = zeros(1,length(R));

for i=1:length(R)
    aR = a./R(i);
    
    
    thr=0.1;
    if gp<=thr
        
        xl=gp*aR^2;
        %  lowL=R^2*aR^2*pi*gp*(-1/xl-xl/4);
        lowL=R(i)^2*aR^2*pi*(-1/aR^2-gp*xl/4);
        
        xu=gp;
        lowU=R(i)^2*aR^2*pi*(-1-gp*xu/4);
        
        
        s_tissue(i)=lowU-lowL;
        %s_tissue=(R^2+gp^2*a^4/4/R^2-a^2-gp^2*a^2/4)*pi;
        
        %  s_tissue=R^2*(1+gp^2*aR^4/4-aR^2-gp^2*aR^2/4)*pi;
        
    elseif gp*aR^2<=thr
        
        func=@(x) besselj(0,x)./x.^2;
        %s_tissue=(-2.125*p+R^2+0.25*p^2/R^2)*pi+pi*p*integral(func,[0.5,gp]);
        
        xl=gp*aR^2;
        lowL=R(i)^2*aR^2*pi*(-1/aR^2-gp*xl/4);
        
        xu=thr;
        lowU=R(i)^2*aR^2*pi*(-gp/xu-gp*xu/4);
        
        s_tissue0=lowU-lowL;
            
        s_tissue(i)=s_tissue0+R(i)^2*aR^2*pi*gp*integral(func,thr,gp);
        
        
    else
        
        func=@(x) besselj(0,x)./x.^2;  % besselj(0,x) ~ (1-x^2/4) at small x
        s_tissue(i)=R(i)^2*aR^2*pi*gp*integral(func,gp*aR^2,gp);
        
    end
    
    
    
end













