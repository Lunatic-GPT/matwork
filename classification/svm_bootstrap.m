function [accuracy,acc_mean,acc_std] = svm_bootstrap(data,labels,frames,mask,nex,nrep,prefix)
% [accuracy,acc_mean,acc_std] = svm_bootstrap(data,labels,frames,mask,nex,nrep,prefix)
% This file test the accuracy of svm_light classifier using bootstrap
% method. Of all the patterns specified in frames. "nex" patterns will be randomly selected for training and the remaining will be used for testing. 
% The random selection is repeated "nrep" times.    
% data: The name of a 4-D afni data file.  Each subbrik of the data file
% correspond to one spatial pattern.     
% labels:  a text file with labels for the subriks in data.  Each line contains one label.  
%                   The number of labels should equal to the number of
%                   subbriks in data.  labels are 1, -1, or 9999. 
%                   If set as 9999, that subbrik will not be used in either
%                   training nor testing.  
% frames: frames to include in data and labels(one based)
% mask:   mask in afni format.  Voxels with non-zero values in the mask will be used in the patterns. 
% nex: number of examples
% nrep: number of repetitions.
% prefix: prefix for the model and testing result
% accuracy: an 1 by nrep array of testing accuracy.

lb = load(labels);
lb_frames = lb(frames);
frames = frames(lb_frames<999);
nfr = length(frames);

if nex >= nfr
    error('Number of examples too big');
end

fr_train = zeros(nex,nrep);
fr_test = zeros(nfr - nex,nrep);

for i=1:nrep
      ind = randperm(nfr);
      fr_train(:,i) = ind(1:nex);
      fr_test(:,i) = ind(nex+1:end);
end

accuracy = zeros(1,nrep);

for i=1:nrep
    svm_learn3d(data,frames(fr_train(:,i)),mask,labels,sprintf('%s_svm_model_bs%d.1D',prefix,i));
    svm_classify3d(data,frames(fr_test(:,i)),mask,labels,sprintf('%s_svm_model_bs%d.1D',prefix,i),sprintf('%s_svm_out_bs%d.1D',prefix,i));
    accuracy(i) = get_accuracy(sprintf('%s_svm_out_bs%d.1D',prefix,i),lb(frames(fr_test(:,i))));
end

acc_mean = mean(accuracy);
acc_std = std(accuracy);

function ind_c=ind_remain(n,ind_arr)
       count = 1;
       ind_c = zeros(1,n-length(ind_arr));
       for i=1:n
           if ~any(ind_arr == i)
               ind_c(count) = i;
               count =  count + 1;
           end
       end

function a = get_accuracy(file,label)

val = load(file);
correct = 0;
for i=1:length(val)
    if val(i)*label(i) > 0
        correct = correct+1;
    end
end
a = correct/length(val);

