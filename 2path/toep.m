function res=toep

d=mfilename('fullpath');
d=fileparts(d);
res=fullfile(d,'..','3TData\Epilepsy3T');

if nargout==0
cd(res);
end

