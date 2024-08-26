function dcm2nii(dpattern,prefix,dout,extract_dcm)
% files will be saved in dout

a=name4pat(dpattern,true);
a=str2cell(a);

for i=1:length(a)
    if ~exist(a{i},'dir')
        continue;
    end
    dname=a{i};
    if ~exist('extract_dcm','var')
        extract_dcm=false;
    end
    
    if ~exist('dout','var')
        dout=pwd;
    end
    if ~exist('prefix','var') || isempty(prefix)
        prefix_new=dname;
    else
        prefix_new=prefix;
    end
    
    if ~exist([filename(prefix_new),'.nii.gz'],'file')
        if isunix
         cmd=sprintf('dcm2niix_afni  -z y -f %s -o "%s" "%s"',filename(prefix_new),dout,dname);   
        else
         cmd=sprintf('dcm2niix  -z y -f %s -o "%s" "%s"',filename(prefix_new),dout,dname);
        end
        
        system(cmd);
        
        if extract_dcm
            extp(dname);
            cdp(dname);
        end
        
    end
    
end