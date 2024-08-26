function dim=dimxform(dim,xform,forw_back)
% forw_back: 0 forward; 1 backward; 
% 0 forward means applying the xform in
% sequence to d.
% backward means apply the inverse tranform to the already xformed data to
% get the original data
% xform: cell array of transforms: 'xy','xz','yz','flipx','flipy','flipz';



if forw_back==0
    for i=1:length(xform)
        if strcmp(xform{i},'xy')
            dim=dim([2,1,3]);
        elseif  strcmp(xform{i},'xz')
            dim=dim([3,2,1]);
        elseif    strcmp(xform{i},'yz')
            dim=dim([1,3,2]);
        end
        
    end
    
else
    
    for i=1:length(xform)
        if strcmp(xform{end-i+1},'xy')
            dim=dim([2,1,3]);
        elseif  strcmp(xform{end-i+1},'xz')
            dim=dim([3,2,1]);
        elseif    strcmp(xform{end-i+1},'yz')
            dim=dim([1,3,2]);
        end
        
    end
    
    
end



