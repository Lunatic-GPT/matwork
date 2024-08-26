function out = size_adapt(image)
% size_adapt rescales the input image with size 256*256 to the DCE image
% size 512*512 pixels.


out = zeros(512,512);
% enlarge DTI image
if size(image,1)==256
    for i = 1:256
        for j = 1: 256
            out(2*i,2*j) = image(i,j);
            out(2*i-1,2*j) = image(i,j);
            out(2*i-1,2*j-1) = image(i,j);
            out(2*i,2*j-1) = image(i,j);
        end
    end
else
    disp('Information: Function size_adapt did not perform any transformation');
    out = image;
end