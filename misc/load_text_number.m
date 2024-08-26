function res=load_text_number(fname)
% load a text file; if the file does not exist, then return []
% if the text file is a matrix then output is a matrix;
% if fname is a string of numbers, then return num2str(fname);
% if fname is [], return [];
res=[];
if ~isempty(fname)
    if isa(fname,'char') 
        if exist(fname,'file')       
          res=load(fname);
        else
            
           res=str2num(fname); 
        end
    else
        
        res=fname;
    end
    
end