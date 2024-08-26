function [h,im_matched, theta,x,y]=image_registr_MI(image1, image2, angle, step)
% function [h,im_matched, theta,I,J]=image_registr_MI(image1, image2, angle, step);
%
% Function for registering image 1 and image 2 using Mutual Information
% Image processing toolbox is required for functions IMROTATE and IMCROP
% For users without IP toolbox download file im_reg_mi.m without
% cropping option and with different rotation function
%
% Mutual Information files are taken from:
% http://www.flash.net/~strider2/matlab.htm
%
% Rigid registration - only translation and rotation are corrected
% For each angle of rotation set of translation is checked
% Can use only for translation by setting angle=0
%
% angle - vector of angles to check, for example 
% angle=[-30:2:30] or 
% angle=15;
%
% step - how many pixels to shift in x and y directions for translation check
%
% crop = 0 to eliminate cropping image
% crop=1 to crop the image and save computational time 
% You'll be asked to crop an area out of the original IMAGE2 
% If you have the kowledge of where approximately the matching area is the
% cropping allows you to limit the search to this area and to save
% calculation time.  Otherwise, select crop=0 
%
% OUTPUT:
% im_matched - matched part of image 2
% h - MI for the best theta
% theta - best angle of rotation
% x and y - shifts in the horizontal and vertical directions
%
% MIGHT BE REALLY SLOW FOR A LOT OF ANGLES AND SMALL STEPS!
% TOO LARGE STEPS CAN CAUSE FALSE REGISTRATION
% I RUN IT OVERNIGHT WITH SMALL STEPS TO GET THE BEST POSSIBLE MATCH
% 
% example:
%
% [h,im_matched, theta,I,J]=image_registr_MI(image1, image2, [-10:5:10], 5, 0);
%
% just translation:
% [h,im_matched, theta,I,J]=image_registr_MI(image1, image2, 0, 1, 1);
%
% written by K.Artyushkova
% 10_2003
% 5-6-2009 XZ, modified to register images of the same size.



if isa(image1,'uint16')
    image1=double(image1)/65535*255;
elseif isa(image1,'uint8')
    image1=double(image1);
else 
    error('Wrong image data format.');
end

image2_o = double(image2);
if isa(image2,'uint16')
    image2=double(image2)/65535*255;
elseif isa(image2,'uint8')
    image2=double(image2);
else 
    error('Wrong image data format.');
end

[m,n]=size(image1);
[p,q]=size(image2); 

if (m~=p || n~=q)
    error('The input image sizes differ');
end

b=length(angle);

%method = 'Normalized';
method = 'Standard';
h = zeros(b,2*step+1,2*step+1);
for k=1:b
    J = imrotate(image2, angle(k),'bilinear','crop'); %rotated cropped IMAGE2
    image21=round(J);
    
    for i=-step:step
        for j=-step:step
            
            if(i<0)  % image 2 shift to left
                xl1 = 1;
                xl2 = -i+1;
                xu1 = n +i;
                xu2 = n;
            else % image 2 shift to right
                xl1 = i+1;
                xl2 = 1;
                xu1 = n;
                xu2 = n-i;      
            end
            
            if(j<0)  % image 2 shift up
                yl1 = 1;
                yl2 = -j+1;
                yu1 = m +j;
                yu2 = m;
            else % image 2 shift down
                yl1 = j+1;
                yl2 = 1;
                yu1 = m;
                yu2 = m-j;      
            end
                im1=image1(yl1:yu1,xl1:xu1); 
                im2=image21(yl2:yu2,xl2:xu2); 
                 
                im2=round(im2); 
                im1=round(im1);
                h(k,step+1+i,step+1+j)=mi2(im1,im2,method); % calculating MI
        end
    end
end
  

[a, b] = max(h(:));% finding the max of MI and indecises
[K,I,J] = ind2sub(size(h),b);

y = I-step-1;
x = J-step-1;

theta=angle(K);
im_rot = imrotate(image2_o, theta,'bilinear','crop');

im_matched=circshift(im_rot,[y,x]);





function h=mi2(image_1,image_2,method)
% function h=MI2(image_1,image_2,method)
%
% Takes a pair of images and returns the mutual information Ixy using joint entropy function JOINT_H.m
% 
% written by http://www.flash.net/~strider2/matlab.htm


a=joint_h(image_1,image_2); % calculating joint histogram for two images
[r,c] = size(a);
b= a./(r*c); % normalized joint histogram
y_marg=sum(b,1); %sum of the rows of normalized joint histogram
x_marg=sum(b,2);%sum of columns of normalized joint histogran

Hy=0;
for i=1:c;    %  col
      if( y_marg(i)==0 )
         %do nothing
      else
         Hy = Hy + -(y_marg(i)*(log2(y_marg(i)))); %marginal entropy for image 1
      end
end
   
Hx=0;
for i=1:r;    %rows
      if( x_marg(i)==0 )
         %do nothing
      else
         Hx = Hx + -(x_marg(i)*(log2(x_marg(i)))); %marginal entropy for image 2
      end   
end
h_xy = -sum(sum(b.*(log2(b+(b==0))))); % joint entropy

if strcmp(method,'Normalized')
h = (Hx + Hy)/h_xy;% Mutual information
else
h = Hx + Hy - h_xy;% Mutual information
end

function h=joint_h(image_1,image_2)
% function h=joint_h(image_1,image_2)
%
% takes a pair of images of equal size and returns the 2d joint histogram.
% used for MI calculation
% 
% written by http://www.flash.net/~strider2/matlab.htm


rows=size(image_1,1);
cols=size(image_1,2);
N=256;

h=zeros(N,N);

for i=1:rows;    %  col 
  for j=1:cols;   %   rows
    h(image_1(i,j)+1,image_2(i,j)+1)= h(image_1(i,j)+1,image_2(i,j)+1)+1;
  end
end


