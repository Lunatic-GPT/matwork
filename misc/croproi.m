function data=croproi(data,marg)
% [l,r,u,d]
    
if ~exist('marg','var')
    marg=0;
end

    thr=1;
    i=1;
    while ~any(data(:,i,:)<thr)    
        i=i+1;
    end
    data(:,1:i-1-marg,:) = [];

    i=0;
    while ~any(data(:,end-i,:)<thr)
        i=i+1;
    end    
    data(:,end-i+1+marg:end,:)=[];


    i=1;
    while ~any(data(i,:,:)<thr)    
        i=i+1;
    end
    data(1:i-1-marg,:,:) = [];

    i=0;
    while ~any(data(end-i,:,:)<thr)
        i=i+1;
    end
    
    data(end-i+1:end+marg,:,:)=[];
    