function wash_in

% Wash_in calculates the relative change between the first two DCE
% images of our DCE series. 
% Finally we can display that way the wash-in process.

load 'C:\Documents and Settings\hahntobi\My Documents\MATLAB\DICOM files\image1';
image1 = image;
load 'C:\Documents and Settings\hahntobi\My Documents\MATLAB\DICOM files\image2';
%image2 = image2;
mean1 = mean(mean(image1));
image1(image1<0.1*mean1)=0;
quot = zeros(size(image1,1),size(image1,2));
quot(image1~=0) = (image2(image1~=0)-image1(image1~=0))./image1(image1~=0);
figure; 
imshow(quot)
colormap(jet)
caxis([0 2])
axis([139,176,267,304])