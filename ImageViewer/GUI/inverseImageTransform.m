function d=inverseImageTransform(d,handles)


it=getappdata(handles.figure1,'ImageTransforms');

for i=1:length(it)
    
    tag=it(end-i+1);
    if strcmp(tag,'flipx_image')        
        d=flip(d,1);
    elseif strcmp(tag,'flipy_image')
        d=flip(d,2);
    elseif strcmp(tag,'flipz_image')
        d=flip(d,3);
    elseif strcmp(tag,'swapxy')
        new_order=[2,1,3,4];
        d=permute(d,new_order);
    elseif strcmp(tag,'swapyz')
        new_order=[1,3,2,4];
        d=permute(d,new_order);
    elseif strcmp(tag,'swapxz')
        new_order=[3,2,1,4];
        d=permute(d,new_order);
    else
        error('unknown transform');
    end
    
end



