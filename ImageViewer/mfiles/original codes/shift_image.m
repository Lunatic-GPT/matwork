function out = shift_image(image,shift, direc)
% This function compensates overall image shifts
% - it works only on 512*512 images
% - shift is the number of pixels you want the image to be shift
% - direc determines the direction of the shift:
% +1 positive vertical direction (Matlab convention)
% -1 negative vertical direction
% +2 positive horizontal direction
% -2 negative horizontal direction
% usage: % image =
% shift_image(shift_image(image,shift_y_abs,direc_y),shift_x_abs,direc_x);

if size(image,1) ~= 512
    disp('ERROR! Shift_image only works on 512*512 images!');
end

out = image;
out2 = zeros(512,512);

if direc == 1
    for k = 1:shift
        for i = 1:512
            for j=1:512
                if i~=1
                    out2(i,j)=out(i-1,j);
                else
                    out2(i,j)=out(512,j);
                end
            end
        end 
        out = out2;
    end
elseif direc == -1
    for k = 1:shift
        for i = 1:512
            for j=1:512
                if i~=512
                    out2(i,j)=out(i+1,j);
                else
                    out2(i,j)=out(1,j);
                end
            end
        end 
        out = out2;
    end
elseif direc == 2
   for k = 1:shift
        for i = 1:512
            for j=1:512
                if j~=1
                    out2(i,j)=out(i,j-1);
                else
                    out2(i,j)=out(i,512);
                end
            end
        end 
        out = out2;
   end
elseif direc == -2
   for k = 1:shift
        for i = 1:512
            for j=1:512
                if j~=512
                    out2(i,j)=out(i,j+1);
                else
                    out2(i,j)=out(i,1);
                end
            end
        end 
        out = out2;
    end
end


        