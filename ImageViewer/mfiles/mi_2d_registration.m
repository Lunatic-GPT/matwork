%% 2-D rigid body Image registration using mutual information. 
% Rigid body 2-D image co-registration (translation and rotation) is performed
% using maximization of mutual information.  
%
% Function is implemented in c-code and compiled using the Matlab compiler
% to minimize computation time.  
%
% Mutual information joint entropy matrix is computed using the Hanning
% windowed sinc function as the kernel of interpolation, which is the HPV
% estimation method [1].  
%
% Maximization of the joint entropy matrix is carried out using Powel's
% Direction set method [2] with original c-code very slightly modified from
% J.P. Moreau's code[3].  
%
% [1] Lu X, et al., Mutual information-based multimodal image registration
% using a novel joint histogram estimation, Comput Med Imaging Graph
% (2008), doi: 10.1016/j.compmedimag.2007.12.001
%
% [2] Numerical Recipes, The Art of Scientific Computing By W.H. Press, 
% B.P. Flannery, S.A. Teukolsky and W.T. Vetterling, Cambridge University 
% Press, 1986
%
% [3] http://pagesperso-orange.fr/jean-pierre.moreau/Cplus/tpowell_cpp.txt
%
% See also: 

%% Syntax
%   [parameters xy xy_0] = mi_hpv_2d(Reference Image,Floating Image)
%
% Input:
%
%   Reference Image: image that will be compared too.  Must be uint8.  Take
%   care to scale image properly.
%
%   Floating Image: image that will be rotated and translated.  Must be
%   uint8.  Take care to scale image properly.
%
% Output:
%
%   parameters: 3x1 Array with the form [DeltaX  DeltaY  DeltaTheta].  
%   Theta is counterclockwise in-plane rotation in radians.  DeltaX/Y
%   are translations in pixels.
%   
%   xy: Optional output.  8x1 Array with the x and y coordinates of the
%   corners of the output matix.
%
%   xy_0: Optional output.  8x1 Array with the x and y coordinates of the
%   corners of the input matrix.
%
%   The output provides the necessary tools for you to translate the Float
%   image, but does not move it for you.  You have the option of using
%   parameters and moving the reference image (using circshift and imrotate
%   for example) or you can use the xy and xy_0 output to more accurately
%   transform the reference image, but its slightly more complicated.  The
%   inclusion of both more has to do with me being unsatisfied with the
%   lack of precision in circshift and imrotate.  I recommend the latter
%   approach.

%% Example
% We will use the Shepp-Logan phantom (phantom.m).  Note we first scale the
% image to have signal intensities in the range 0-255 and then convert to
% uint8.  
function mi_2d_registration()
Ref = imread('01.tiff'); % same as: Ref = uint8(phantom.*255);
Float = imread('02.tiff');
data_dir = '../Data';
dir_struct = dir(fullfile(data_dir, '*X1'));
    
    cut_off_rect_pos = fullfile(data_dir,dir_struct.name, 'cut_off_rect_pos.mat');
    load(cut_off_rect_pos, 'rect_pos');
        x1 = round(rect_pos(1));
        x2 = round(rect_pos(3));
        y1 = round(rect_pos(2));
        y2 = round(rect_pos(4));
        
Ref_cp = Ref(y1:(y1+y2),x1:(x1+x2));
Float_cp = Float(y1:(y1+y2),x1:(x1+x2));
     
motion_check(Ref_cp,Float_cp);

[h,im_matched, theta,I,J]=image_registr_MI(Ref_cp, Float_cp,0,,crop);



