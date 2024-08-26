function mat2afni(d,afni,prefix)

if isa(d,'char')
    
    
    if ~exist('prefix','var')
        prefix=strtok(d,'.');
    end
    
    d=ri(d);
    
end

if ~exist('prefix','var')
    prefix=strtok(d,'.');
end

[err,info]=BrikInfo(afni);

WriteBrikEZ(d,info,'',prefix,'');



