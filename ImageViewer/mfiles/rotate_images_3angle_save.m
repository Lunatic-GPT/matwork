function d2=rotate_images_3angle_save(d,ph,prefix)
% d2=rotate_images(d,ph,preifx)
% ph: 3 angles in degrees
% prefix: output file name.  If not exist, then will take from d (d has to
% be a file name).
% 

if isa(d,'char')
    if ~exist('prefix','var')
        
    prefix=strtok2(d,'.');
    prefix=strtok(prefix,'+');
    
    ph_str='';
    for i=1:length(ph)
        ph_str=[ph_str,'_',num2str(ph(i))];        
    end
    
    prefix=[prefix,ph_str,'.mat'];
    end
    
     
     d=ri(d);
end

d2=d;

for i=1:3
    d2=rotate_images(d2,i,ph(i));
end


save(prefix,'d2');