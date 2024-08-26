if exist('gb_unique','var') || exist('dbstack_str','var') 
    error('variable gb_unique already exist');
end
gb_unique=whos;
gb_unique=structarray(gb_unique,'bytes');
gb_unique=sum(gb_unique)/(1024)^3;
dbstack_str=dbstack();
fprintf('%s: memory = %f GB\n',dbstack_str(2).name,gb_unique);
clear gb_unique dbstack_str;

