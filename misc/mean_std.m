function mn_sd=mean_std(d,dim,do_squeeze)

if ~exist('do_squeeze','var')
    do_squeeze=false;
end

if ~exist('dim','var')%assuming a vector
    dim=2;
    d=reshape(d,[1,length(d)]);
end
    
mn=mean(d,dim);
sd=std(d,[],dim);

if do_squeeze
    mn=squeeze(mn);
    sd=squeeze(sd);
end
mn_sd=[mn,sd];