function res=fmt_num2name(fmt,varargin)

str={};
for i=1:length(varargin)
   str{i}=num2str(varargin{i});
    
end
res=sprintf(fmt,str{:});

