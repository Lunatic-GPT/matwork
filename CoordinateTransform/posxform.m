function pos=posxform(pos,xform,sz,forw_back)
% forw_back: 0 forward; 1 backward; 
% 0 forward means applying the xform in
% sequence to d.
% backward means apply the inverse tranform to the already xformed data to
% get the original data
% xform: cell array of transforms: 'xy','xz','yz','flipx','flipy','flipz';
% pos: n*3


if forw_back==0
   i_array=1:length(xform);
else
    i_array=length(xform):-1:1;
end
   

    for i=i_array
        if strcmp(xform{i},'xy')
            pos=pos(:,[2,1,3]);
        elseif  strcmp(xform{i},'xz')
            pos=pos(:,[3,2,1]);
        elseif    strcmp(xform{i},'yz')
            pos=pos(:,[1,3,2]);
        elseif strcmp(xform{i},'flipx')
            pos(:,1)=sz(1)-pos(:,1)+1;
        elseif strcmp(xform{i},'flipy')
            pos(:,2)=sz(2)-pos(:,2)+1;
        elseif strcmp(xform{i},'flipz')
            pos(:,3)=sz(3)-pos(:,3)+1;
        end
        
        
        
    end
    
    

