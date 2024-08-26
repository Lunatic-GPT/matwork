function [data,rows,cols]=crop_zero(data)

sz=size(data);
    i=1;
    thr=1;
    while ~any(data(:,i,:)>0)    
        i=i+1;
    end
    data(:,1:i-1,:) = [];

    cols(1)=i;
    
    i=0;
    while ~any(data(:,end-i,:)>0)
        i=i+1;
    end    
    data(:,end-i+1:end,:)=[];
   cols(2)=sz(2)-i;

    i=1;
    while ~any(data(i,:,:)>0)    
        i=i+1;
    end
    data(1:i-1,:,:) = [];

    rows(1)=i;
    i=0;
    while ~any(data(end-i,:,:)>0)
        i=i+1;
    end
    
    data(end-i+1:end,:,:)=[];
    
    rows(2)=sz(1)-i;
    
    cols=cols(1):cols(2);
    rows=rows(1):rows(2);
    
    
    
    