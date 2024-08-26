function val = parValArray(fid_prefix,name)
%parValArray(fid_prefix,name)

array = readPar([fid_prefix,'.fid'],'array');
[name1,name2]=strtok(array(2:end-1),',');


if isempty(name2) 
    if strcmp(name1,name)
        val = readPar(fid_prefix,name);
    else
        error('Parameter %s was not arrayed',name);
    end
    
else
    name2 = name2(2:end);
    v1 = readPar(fid_prefix,name1);
    v2 = readPar(fid_prefix,name2);
        
    if strcmp(name1,name)
      x = meshgrid(v1,v2);
      val = x(:);
    elseif strcmp(name2,name)
      x = meshgrid(v2,v1);
      x = x';
      val = x(:);  
    else
        error('Parameter %s was not arrayed',name);
    end
end


    
      
     