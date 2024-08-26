function vmask_Frangi2D(params)

m=get_fpattern(params,'brain mask')
mask=ri('mask_wm_MID35.mat');
 
 
 
 %%
% .FrangiScaleRange : The range of sigmas used, default [1 8]
%       .FrangiScaleRatio : Step size between sigmas, default 2
%       .FrangiBetaOne : Frangi correction constant, default 0.5
%       .FrangiBetaTwo : Frangi correction constant, default 15
%       .BlackWhite : Detect black ridges (default) set to true, for
%                       white ridges set to false.
%       .verbose : Show debug information, default true
c=ri('Phase_MID35_mean_detrend.mat');

options.FrangiScaleRange= [0.5,2];
options.FrangiScaleRatio = 0.5;
options.BlackWhite=false;
d=FrangiFilter2D(c,options,[],mask);
d2=d;
d2(mask==0)=0;
save Frangi d2
%%









