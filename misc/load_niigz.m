function d=load_niigz(name)


        [nii_name,suf]=strtok2(name,'.');
        
        if strcmp(suf,'.gz')
        gunzip(name);       
        else
            unzip(name);       
        end
        d=load_nii(nii_name);
        
        delete(nii_name);
        
        
