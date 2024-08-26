function [img1_new,img2_new]=select_pixels4ScatterPlot(img1,img2)

% select pixels between 2 lines in a scatter plot of img1 and img2
% the line should be drawn with gline; do not use insert which cannot be
% found by findobj.
h=findobj(gcf,'Type','Line');

x=zeros(2,2); % 2 points*2 lines
y=zeros(2,2);  % 

c=1;
for i=1:length(h)
   
    
    xdata=get(h(i),'XData');
    
    
    ydata=get(h(i),'YData');
    
    if length(xdata)==2   
        x(:,c)=xdata;
        y(:,c)=ydata;
        c=c+1;
    end
end

if c~=3
    error('Number of lines not equal to 2');
end

    y1=get_liney(x(:,1),y(:,1),img1);
    y2=get_liney(x(:,2),y(:,2),img1);
    
    m=(img2>=y1 & img2<=y2) | (img2>=y2 & img2<=y1);
    
    img1_new=img1;
    img1_new(~m)=0;
    
    img2_new=img2;
    img2_new(~m)=0;
    
    
    


function y= get_liney(px,py,x)

% px and py are the x and y coordinates of two points; 1*2


k=diff(py)/diff(px);

b=py(1)-k*px(1);

y=k*x+b;
