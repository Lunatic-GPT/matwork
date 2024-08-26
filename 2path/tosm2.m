function d=tosm2


d=mfilename('fullpath');
d=fileparts(d);
d=fullfile(d,'..','simulation','susceptibility_PartialVolume');
if nargout==0
  cd(d);   
end