function do_dcm2nii(substring)
% do_dcm2nii(substring)
% after run can use move_nii.m to move all the nii files into one folder

a=dir2('*');
           
for i=1:length(a)
    if a(i).isdir
        if ~isempty(strfind(a(i).name,substring))
            dcm2niix(a(i).name);     
        else
            cd(a(i).name);
            if exist('dest','var')
            else
             do_dcm2nii(substring);
            end
            cd('..');
        end
    
    end
    
   
    
end

function dcm2niix(dname)
 cmd=sprintf('dcm2niix -o %s -f %s %s ',pwd,dname,dname);
 unix(cmd);
 
 
 
 
          
