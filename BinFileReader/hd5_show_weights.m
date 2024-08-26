function [w,name]=hd5_show_weights(fname)

in=h5info(fname);
%     %% optimizer_weights
%     for i=1:41
%         disp(i);
%         in.GroupHierarchy.Groups(2).Datasets(i)
%
%     end
w={};
name={};
%% model_weights
for i=1:length(in.Groups(1).Groups)
    for j=1:length(in.Groups(1).Groups(i).Attributes.Value)
        fprintf('%d - %s/%s \n',i, in.Groups(1).Groups(i).Name,in.Groups(1).Groups(i).Attributes.Value{j});
        name{i,j}=sprintf('%s/%s', in.Groups(1).Groups(i).Name,in.Groups(1).Groups(i).Attributes.Value{j});
        w{i,j}=h5read(fname,name{i,j});
    end
    if isempty(in.Groups(1).Groups(i).Attributes.Value)
        fprintf('%d - %s\n',i, in.Groups(1).Groups(i).Name);
    end
    
    
end

%%

