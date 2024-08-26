function [mn,sd]=mstd(d,dim)

if exist('dim','var')
    mn=mean(d,dim);
    sd=std(d,[],dim);
else
    mn=mean(d);
    sd=std(d);
end

if nargout==1 && length(mn)==1
    mn=[mn,sd];
end
