function applyShift(fname,newdir,fparameter)

if ~exist('fparameter','var')
    fparameter='overlayShiftParameters.mat';
end

d=ri(fname);

tmp=load(fparameter);

name=fieldnames(tmp);

for i=1:length(name)
    
    if strcmp(name{i},'step_dim2')
      d=circshift(d,tmp.step_dim2,2);
        
    elseif strcmp(name{i},'step_dim1')
        d=circshift(d,tmp.step_dim1,1);
    else
        error('unknown field');
    end
    
end

save(fullfile(newdir,fname),'d');

