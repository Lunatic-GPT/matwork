function mat_chvname(mat,varargin)
% mat_chvname(mat,oldname,newname)
tmp=load(mat);

for i=1:nargin/2
eval(sprintf('tmp.%s=tmp.%s;',varargin{(i-1)*2+2},varargin{(i-1)*2+1}));
   
tmp=rmfield(tmp,varargin{(i-1)*2+1});
 
end
save(mat,'-struct','tmp');