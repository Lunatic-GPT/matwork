function res=waitforkey(keys)

while 1
    keydown = waitforbuttonpress;
    if keydown==1
        key=get(gcf, 'CurrentCharacter');
        
        if any(key==keys)
            res=find(key==keys);
            break;
        end
    end
    
    
end
    
       