   
function b1=B1_th2th(y,th)

r=y(2)/y(1);
x=(r+sqrt(r^2+8))/4;
            
b1=acos(x)*90/th*180/pi;
            
            