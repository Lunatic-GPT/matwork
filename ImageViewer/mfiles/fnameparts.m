function [s,t,e]=fnameparts(fname)

fname = strtok(fname,'.');

 [st1,r1]= strtok(fname,{'_'});
 
 s = str2double(st1);
 
 if  isnan(s)
     error([mfilename,': Unknown file name pattern']);
 end
 
 if isempty(r1)
     error([mfilename,': Unknown file name pattern']);
 end
 
 
 [st2,r2]=strtok(r1,{'_'});
 
 if isnan(str2double(st2))
     t = [];
     e = r1(2:end);
   
 else
     
     t = str2double(st2);
     
     if ~isempty(r2)
      
      e = r2(2:end);
     else
         e = [];
     end
 
   
 end
 
 
 
 
 