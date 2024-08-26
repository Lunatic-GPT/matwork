function dname=topvs(id)

d=mfilename('fullpath');
d=fileparts(d);

if ~exist('id','var')
dname=fullfile(d,'..','7TData/PVS_R21');
else

    if isa(id,'char')
        id=str2num(id);
    end
   dname= fullfile(d,'..','7TData/PVS_R21',sprintf('PVS%02d',id));
end
if nargout==0
cd(dname); 
end



