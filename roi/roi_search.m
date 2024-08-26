function rois=roi_search(label)
         
tilde = find(label=='~');

ntd = length(tilde);

if (mod(ntd-3,3)~=0) 
    error('Wrong roi label format');
end
nrois = (ntd-3)/3;

fprintf('%d rois found\n',nrois);
rois = cell(1,nrois);
for i=1:nrois

    len = 0;
    while label(tilde(3+i)-len-1)~=' ' ...
       && label(tilde(3+i)-len-1)~='/' ...
       && label(tilde(3+i)-len-1)~='\'
      len=len+1;
    end
    
    rois{i}= label(tilde(3+i)-len:tilde(3+i)-1);
    fprintf('%s  ',rois{i});
end



fprintf('\n');