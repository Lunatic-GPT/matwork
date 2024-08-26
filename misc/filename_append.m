function res=filename_append(fname,append,before_first_dot)

if ~exist('before_first_dot','var')
  before_first_dot=false;
end

 [d,f,a]=fileparts(fname);

if before_first_dot

 [f1,f2]=strtok(f,'.');
 res=fullfile(d,[f1,append,f2,a]);
 
else
 res=fullfile(d,[f,append,a]);
end