function d=load_untouch_niigz(name)

if ~isa(name,'char')
    d=name;
    return;
end

if ~strcmp(name(end-2:end),'.gz')
    d=load_untouch_nii(name);
    return;
end

[nii_name,suf]=strtok2(name,'.');
prefix=strtok2(nii_name,'.');

if ~exist([prefix,'.nii'],'file')
    if strcmp(suf,'.gz') && exist(name,'file')
        gunzip(name);
    elseif strcmp(suf,'.zip') && exist(name,'file')
        unzip(name);
    end

    d=load_untouch_nii(nii_name);
        

  %  save([prefix,'.mat'],'d');
    pause(0.5);%avoid permission denied error for delete
    delete(nii_name);
else  
    
    d=load_untouch_nii(nii_name);
   %{ 
    a=dir( [prefix,'.mat']);
    b=dir(name);
    
    if a.datenum<b.datenum
       delete( [prefix,'.mat']);
       d=load_untouch_niigz(name);
    else   
       d=ri( [prefix,'.mat'],'','','d');
    end
    %}
end

  