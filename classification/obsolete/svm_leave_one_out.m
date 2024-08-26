function [accuracy,nsamples,acc_mean] = svm_leave_one_out(data,labels,frames,mask,ngroups,prefix)
% [accuracy,nsamples,acc_mean] = svm_leave_one_out(data,labels,frames,mask,ngroups,prefix)
% frames will be divided approximately equally into ngroups.  One of the
% group will be used for test and the rest for training.
% data: file name for 3d pattern
% labels: labels for each subbrik pattern in data
% frames: frames to include in data and labels(one based)
% ngroups: number of groups to divide
% prefix: prefix for the model and testing result
% accuracy: an 1 by ngroups array of testing accuracy.

lb = load(labels);


lb_frames = lb(frames);

frames = frames(lb_frames<999);

nfr = length(frames);


ind = 1:round(nfr/ngroups):nfr;

ind2 = round(nfr/ngroups):round(nfr/ngroups):nfr;

if ind2(end) <nfr
    ind2(end+1) = nfr;
end

accuracy = zeros(1,ngroups);
nsamples = zeros(1,ngroups);

for i=1:ngroups
    
    fr_train = [];
    for j=1:ngroups
        if j~=i
          fr_train = [fr_train,ind(j):ind2(j)];
        end
    end
    
    svm_learn3d(data,frames(fr_train),mask,labels,sprintf('%s_svm_model_loo%d.1D',prefix,i));
    svm_classify3d(data,frames(ind(i):ind2(i)),mask,labels,sprintf('%s_svm_model_loo%d.1D',prefix,i),sprintf('%s_svm_out_loo%d.1D',prefix,i));
    accuracy(i) = get_accuracy(sprintf('%s_svm_out_loo%d.1D',prefix,i),lb(frames(ind(i):ind2(i))));
    nsamples(i) = length(fr_train);
end


acc_mean = mean(accuracy);

function a = get_accuracy(file,label)

val = load(file);
correct = 0;
for i=1:length(val)
    if val(i)*label(i) > 0
        correct = correct+1;
    end
end
a = correct/length(val);

