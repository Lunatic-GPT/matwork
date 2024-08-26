function y=R1rasym_multioffset(b,x)

% b has [kpH6_1,2,3,kpH7_1,2,3]
%xdata [sample,B1,off,con]
% B1 and off are in herz.
% con in mM;
% output is an array of R1r. 

%sample number 'pH = 6 (25 mM)','pH = 7.0 (25 mM)','pH = 6.0 (50 mM)','pH = 7.0 (50 mM)'});



R1a=[ 0.2504,0.2514,0.2987,0.2477];

R1=R1a(x(:,1))';

%{

if x(1) ==1 || x(1)==3
    k=b(5);
else
    k=b(6);
end
%}

%R2=b(1);

R2=b(x(:,1))';
    

of=[1.2,2.1,2.8]*400*2*pi;

con=x(:,4)*1e-3/111*[2,1,1];

Rex=zeros(size(x,1),1);
Rex2=zeros(size(x,1),1);

for i=1:size(x,1)
    
    
    if x(i,1) ==1 || x(i,1)==3
      k=b(1:3);
    else
      k=b(4:6);
    end
    B1=x(i,2)*2*pi;
    
    offset=x(i,3)*2*pi;

        
    Rex(i)=sum(con(i).*k.*of.^2./((of-offset).^2+B1^2+k.^2));
    Rex2(i)=sum(con(i).*k.*of.^2./((of+offset).^2+B1^2+k.^2));
 
end


sinth2 = x(:,2).^2./(x(:,2).^2+x(:,3).^2);

y=R1.*(1-sinth2)+(Rex+R2).*sinth2;
y=y;
tmp=(Rex-Rex2).*sinth2;
y(x(:,3)>0)=tmp(x(:,3)>0);  

