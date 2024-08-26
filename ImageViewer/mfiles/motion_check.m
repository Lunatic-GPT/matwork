function motion_check(f1,f2)

handle = figure;
colormap(gray);
subplot(1,2,1);

imageProp = { 'ButtonDownFcn'};
    
    imageVal = { 'copy_rect()' };
    
    imagesc(f1);

h_rect = imrect(gca, [1 1 50 50]);

% click the mouse to copy hrect from left.

subplot(1,2,2);
imagesc(f2, imageProp,imageVal); colormap(gray);
    
 

        
        
       