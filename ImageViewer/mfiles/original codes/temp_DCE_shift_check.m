function temp_DCE_shift_check

close all
SID = input('Please denote the subject ID: ','s');
temp = ['\\vali\hahntobi\matlab\patient_images\' SID '\DCE_temp\01.mat'];
load(temp,'image');
DTI = size_adapt(image);
figure; 
imshow(DTI,[]);
if nargin==0, fig = gcf; end
units = get(fig,'units');
set(fig, 'units', 'normalized', 'outerposition', [0 0 1 1]);
set(fig, 'units', units);
temp = ['\\vali\hahntobi\matlab\patient_images\' SID '\DCE_temp\02.mat'];
load(temp,'image')
DCE = image;
figure; 
imshow(DCE,[]);
if nargin==0, fig = gcf; end
units = get(fig,'units');
set(fig, 'units', 'normalized', 'outerposition', [0 0 1 1]);
set(fig, 'units', units);
temp = ['\\vali\hahntobi\matlab\patient_images\' SID '\DCE_temp\03.mat'];
load(temp,'image');
DTI = size_adapt(image);
figure; 
imshow(DTI,[]);
if nargin==0, fig = gcf; end
units = get(fig,'units');
set(fig, 'units', 'normalized', 'outerposition', [0 0 1 1]);
set(fig, 'units', units);