function wash_out

% Wash_out calculates the slopes of each pixel in the time series of image
% 3-6. 
%%%Be careful with using this function from the .m-file. Better copy the
%%%whole file to the command window. There are some strange reasons why the
%%%program doesn't work as a .m-file properly (see below).

load 'C:\Documents and Settings\hahntobi\My Documents\MATLAB\DICOM files\image3';
image3 = image; % Executing this from the .m-file gives wrong results! Execute this code in the command window!
load 'C:\Documents and Settings\hahntobi\My Documents\MATLAB\DICOM files\image4';
image4 = image;
load 'C:\Documents and Settings\hahntobi\My Documents\MATLAB\DICOM files\image5';
image5 = image;
load 'C:\Documents and Settings\hahntobi\My Documents\MATLAB\DICOM files\image6';
image6 = image;

s = zeros(512,512);
my_array = cat(3,image3, image4, image5, image6);
for i = 1: 512
    for j = 1:512
        temp = (squeeze(my_array(i,j,:)))';
        fit_vector=polyfit([1 2 3 4],temp,1);
        slope = fit_vector(1);
        s(i,j) = slope;
    end
end
save('C:\Documents and Settings\hahntobi\My Documents\MATLAB\temporary saved variables\wash_out', 's');

s_norm = s/180;
s_rad = atan(s_norm)*(2*90/pi); figure; imshow(s_rad,[])
caxis([-90 90])
jet_inv = jet;
jet_inv = flipud(jet_inv);
colormap(jet_inv)
axis([139,176,267,304])