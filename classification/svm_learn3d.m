function svm_learn3d(fname,frames,fmask,category,model)
%svm_learn3d(fname,frames,fmask,category,model)
% input fields: 
% fname: The name of a 4-D afni data file.  Each subbrik of the data file
% correspond to one spatial pattern.     
% frames:  indices of subbriks to be used for training.  1 based.
% fmask:             mask in afni format.  Voxels with non-zero values in the mask will be used in the patterns. 
% category:  a text file with labels for the subriks in data.  Each line contains one label.  
%                   The number of labels should equal to the number of
%                   subbriks in data.  labels are 1, -1, or 9999. 
%                   If set as 9999, that subbrik will not be used in training.  
% model:  name of svm model output file.
%  Details on the format of Model and test output files can be found on http://www.cs.cornell.edu/People/tj/svm_light/

tic;
svm_afni2txt(fname,fmask,category,'svm_learn3d_temp.1D',frames);

disp('Training ... ');

labels = load(category);
if any(labels<0)
 cmd = sprintf('svm_learn svm_learn3d_temp.1D %s',model);
else
 cmd = sprintf('svm_multiclass_learn svm_learn3d_temp.1D %s',model);
end

unix(cmd);

disp([mfilename,' finished in ',num2str(toc),' seconds.']);
