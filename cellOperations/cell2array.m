function [res,mn,id]=cell2array(c,dim,nb,ne)
% [res,mn]=cell2array(c,dim,nb,ne)
% mn gives the mean of each cell element
% nb: the beginning nb elements are skipped
% ne: the last ne elements are skipped

res=[];
id=[];
if ~exist('nb','var')
    nb=0;
end

if ~exist('ne','var')
    ne=0;
end

for i=1:length(c)
    if ~exist('dim','var') || isempty(dim)
        res=[res;vec(c{i}(nb+1:end-ne))];
        mn(i)=mean(c{i}(nb+1:end-ne));
        id=[id;i*ones(length(c{i}(nb+1:end-ne)),1)];
    else
        ctmp=get_part(c{i},dim,nb,ne);
        res=cat(dim,res,ctmp);
        mn(i,:,:,:)=mean(ctmp,dim);

        id=[id;i*ones(size(ctmp,dim),1)];
    end
    
end

mn=squeeze(mn);


function res=get_part(d,dim,nb,ne)

p=1:ndims(d);

p(1)=dim;
p(dim)=1;


d=permute(d,p);

d=d(nb+1:end-ne,:,:,:,:,:);

res=permute(d,p);

