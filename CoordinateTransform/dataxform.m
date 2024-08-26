function d=dataxform(d,xform,forw_back)
% forw_back: 0 forward; 1 backward; 
% 0 forward means applying the xform in
% sequence to d.
% backward means apply the inverse tranform to the already xformed data to
% get the original data
% xform: cell array of transforms: 'xy','xz','yz','flipx','flipy','flipz';



if forw_back==0
    for i=1:length(xform)
        if strcmp(xform{i},'xy')
            d=permute(d,[2,1,3,4]);
        elseif  strcmp(xform{i},'xz')
            d=permute(d,[3,2,1,4]);
        elseif    strcmp(xform{i},'yz')
            d=permute(d,[1,3,2,4]);
        elseif    strcmp(xform{i},'flipx')
            d=flip(d,1);
        elseif    strcmp(xform{i},'flipy')
            d=flip(d,2);
        elseif    strcmp(xform{i},'flipz')
            d=flip(d,3);
        end
        
    end
    
else
    
    for i=1:length(xform)
        if strcmp(xform{end-i+1},'xy')
            d=permute(d,[2,1,3,4]);
        elseif  strcmp(xform{end-i+1},'xz')
            d=permute(d,[3,2,1,4]);
        elseif    strcmp(xform{end-i+1},'yz')
            d=permute(d,[1,3,2,4]);
        elseif    strcmp(xform{end-i+1},'flipx')
            d=flip(d,1);
        elseif    strcmp(xform{end-i+1},'flipy')
            d=flip(d,2);
        elseif    strcmp(xform{end-i+1},'flipz')
            d=flip(d,3);
        end
        
    end
    
end



