function res3=readdPar(dname,par,allfile)
%res=readdPar(dname,par,allfile)

d=dir2(dname);
warning off;
if ~exist('allfile','var')  || ~allfile
    
    for i=1:length(d)
        if isdir(fullfile(dname,d(i).name))
            continue;
        end
        try
        in=dicominfo(fullfile(dname,d(i).name));
        
        break;
        catch
            continue;
        end
    end
    if exist('par','var')
        res3=getfield(in,par);
    else
        
        disp(in);
    end

else
    
    res3={};
    for i=1:length(d)
     %   disp(i);
        in=dicominfo(fullfile(dname,d(i).name));
        try
            res3{i}=getfield(in,par);
        catch
            res3{i}='';
        end
        
    end
    
    
    if isa(res3{1},'char')
       return; 
    end
    try
     res3=cell2mat(res3);
    catch
        
    end
    
end

    
%disp(res3);
