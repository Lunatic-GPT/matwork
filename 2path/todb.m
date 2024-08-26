function d=todb


d='C:\Users\zongx\OneDrive - University of North Carolina at Chapel Hill\Projects_ongoing\Diabetes';
if nargout==0
    cd(d);
end

%{
d=mfilename('fullpath');
d=fileparts(d);
if nargout==0
    cd(d);
    cd('..');
else
    dname=fullfile(d,'..');
end
%}