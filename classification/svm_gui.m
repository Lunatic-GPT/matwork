function svm_gui
% GUI for two-way support vector machine analysis using afni dataset.
% input fields: 
% data: The name of a 4-D afni data file.  Each subbrik of the data file
% correspond to one spatial pattern.     
% category labels:  a text file with labels for the subriks in data.  Each line contains one label.  
%                   The number of labels should equal to the number of
%                   subbriks in data.  labels are 1, -1, or 9999. 
%                   If set as 9999, that subbrik will not be used in either training nor testing.  
% frames (training):  indices of subbriks to be used for training.  1 based.
% frames (testing):  indices of subbriks to be used for testing.  1 based.
% mask:             mask in afni format.  Voxels with non-zero values in the mask will be used in the patterns. 
% Model:  name of svm model output file.
% Test output: name of the test result file.
%  Details on the format of Model and test output files can be found on http://www.cs.cornell.edu/People/tj/svm_light/

p = parameter('svm_light (Support Vector Machine) interface');

p = add(p,'filename','data','');
p = add(p,'filename','category labels','');

p = add(p,'int','frames(training)',[]);
p = add(p,'int','frames(test)',[]);

p = add(p,'filename','mask','');
p = add(p,'filename','model','svm_model.1D');

p = add(p,'filename','Test Output','svm_classify_output.1D'); 
p = add(p,'pop-up menu','action',{'train+test','train','test'});
p = add(p,'button','Go','svm_gui_callback(params)');
p = add(p,'button','Close','close');
p = add(p,'button','Help','help svm_gui');
p = parametergui(p);

