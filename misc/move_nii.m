function move_nii(dpattern,nii_dir,prefix)
% the nii files will be named

a=dir(dpattern);

for i=1:length(a)
    
    if ~a(i).isdir
        continue;
    end
    
    b=dir(fullfile(a(i).name,'*.nii'));
    
    for j=1:length(b)
        
        movefile(fullfile(b(j).folder,b(j).name),fullfile(nii_dir,sprintf('%s_%s.nii',prefix,a(i).name)));
    end
end
    