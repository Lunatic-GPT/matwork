function BLAT
% -----------------------------------------------------------
% Copyright Department of Radiology
% Michigan State University, 2007/2008
%
% BLAT is a software designed to analyse breast lesion
% It is written and tested for UNIX systems.
%
%
% ------------------
% -- INSTALLATION --
% ------------------
% In order for BLAT to work, you simply need to copy the folders "Easy GUI", "Image Gallery" 
% and "m-files" into a directory of your choice. Add it with all subdirectories to your
% Matlab path.
% 
%
% ------------------------
% -- GENERAL INFORMATION--
% ------------------------
% BLAT allows you to create images from dicom data, find breast lesion
% margins and calculate specific values on the created lesion, tissue 1 and
% tissue 2 areas. BLAT is split into two parts, a combined DCE/DTI analysis
% that is based on max. 3 slices and a detailed DCE analysis, based on
% multiple slices. 
%
%
% -----------
% -- USAGE -- 
% -----------
%
% DOWNLOAD IMAGES
% ---------------
% Press 'DOWNLOAD IMAGES' to create the patient images of interest from the
% DICOM files available on the central server. 
% Adjust source and destination folders if necessary and make sure to use
% the right slice numbers.
% 'DCE iteration' denote the total number of DCE slices in one 3d volume image.
%
% 
% SHOW IMAGES/START ANALYZE
% -------------------------
% Press 'SHOW IMAGES/Define ROIs' to open the imageGallery GUI. 
% Press 'Select images' to select images to display or press 'Select all slices
% at t=2' to display all the slices in a directory.  In the latter case, only the first 
% post-contrast images will be displayed.  
% Before selecting images to display, choose the type of scaling you want.
% Auto constant: all the images will use the same colormap and white correspons 
%                to the maximum value of all the images.
% Auto individual: Different colormaps will be used for the images and
%                 white in each image corresponds to the maximum value of that image.
% Manual constant: All the images will use the same colormap and you can
%                 sepcify the image value for white manually.
% After the images are displayed, press 'Add crops' button to add a small square to all 
% the images.  Adjust the position and size of the square in one of the images.  Then click
% 'Copy crops' to copy the new crop square to the other images.  Click
% 'Clear crops' to remove the squares.  
% Press 'Load ROI' to load a predefined roi and inner and outer boundaried to all the images.  
% Press 'Clear ROI' to clear the rois.
% 
% Left click on an image outside the crop region to display the image
% in a new window where you can define rois and save rois.  

% During the next steps, you will draw the rough lesion area and the inner
% and outer boundaries. Then BLAT will detect the lesion margins.
% 
% 3D ROI
% ----------------------------------------------------------------------
%  the ROIs in the prevoius step are detected slice by slice.  Here you
%  can use the auto detection tool to detect a 3-dimensional roi, using the 
% above rois as a starting point.  You will be able to compare these two 
% different rois and choose which one to use in 'Detailed DCE study'

% SHIFT DOWNLOADED DCE IMAGES
% ---------------------------
% By all means you have to perform this step, EVEN IF THERE IS NO SHIFT IN THE DCE SERIES!
% (Not only the images will be shifted, but also kinetic images like the
% wash-out or the wash-in image will be created.)
% The first after contrast image will be used as the base image.  
% The motion correction parameters will be calculated automatically.   
% The data directory is the 'Data' directory created in the previous step.
%
%
% 
%  The same is valid for the
% "Detailed DCE study", where you might have downloaded much more slices
% than you actually analyze. The files for the individual slices are
% number, however, from 1 to the maximum number. The "first slice to be
% analyzed" denotes the slice with the specific number that you can also
% find by looking at the image names in your image directory.
% 
%
% DETAILED DCE STUDY
% ------------------
% use 3dROI: whether to use the roi detected in 3D (from above '3D ROI' step)
% "Right cut-off" and "Left cut-off" denote the two slope values that separate 
% the lesion into three regions (wash-out, plateau, persistent). The values
% are to be given in degree.


BLAT=parameter('BLAT v.1.0, MSU');

BLAT=add(BLAT,'button','Download images','download_images');
BLAT=add(BLAT,'button','Show images/Define ROIs','imageGallery');

%BLAT=add(BLAT,'button','3D ROI','roiAdjust3D');
%BLAT=add(BLAT,'button','Motion correction/Calculate wash in&out', 'auto_shift_DCE_images');
BLAT=add(BLAT,'button','Calculate wash in&out', 'Analyze');
BLAT=add(BLAT,'button','Show wash in&out', 'show_wash_in_out');
%BLAT=add(BLAT,'button','Detailed DCE study','combined_analysis_DCE_study_prefunction');
BLAT = add(BLAT,'button','Quit','close');

BLAT=parametergui(BLAT);

