function svm_gui_callback(p)


f = get(p,'data');
labels = get(p,'category labels');
fr_train = get(p,'frames(training)');
fr_test = get(p,'frames(test)');

mask = get(p,'mask');
model = get(p,'model');
f_out = get(p,'Test Output');

action = get(p,'action');

if strcmp(action,'train')
    svm_learn3d(f,fr_train,mask,labels,model);
elseif strcmp(action,'test')
    svm_classify3d(f,fr_test,mask,labels,model,f_out);
else
    svm_learn3d(f,fr_train,mask,labels,model);
    svm_classify3d(f,fr_test,mask,labels,model,f_out);
end




