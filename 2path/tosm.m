function d=tosm


d=mfilename('fullpath');
d=fileparts(d);
d=fullfile(d,'..','simulation');
if nargout==0
  cd(d);   
end