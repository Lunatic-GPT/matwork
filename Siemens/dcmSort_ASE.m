function [name_sort,ns]=dcmSort_ASE(dname, isold)

% old_new: if old is true; the images before FFT along z-direction was not
% save.

extp(dname);
ns=readsPar([dname,'.pro'],'sSliceArray.lSize');

com=readdPar(dname,'ImageComments',true);

dir_str=dir(dname);

dnum=[];
name={};
for i=1:length(dir_str)
    
    num=str2num(dir_str(i).name);
    if ~isnan(num)
        dnum(end+1)=num;
        name{end+1}=dir_str(i).name;
    end
end


if isempty(dnum)
 
 d=ri(dname,1);
 name_sort=[];
else
    
    [tmp,ind]=sort(dnum);
    
     name_sort=name(ind);

    for i=1:length(dnum)
     
        if i==1
            d=dicomread(fullfile(dname,name{ind(i)}));
        else
            d2=dicomread(fullfile(dname,name{ind(i)}));
            d=cat(3,d,d2);
        end  
    end
    

end

d=squeeze(d);
sz=size(d);

%%
nte=readsPar([dname,'.pro'],'alFree[10]');
te=zeros(1,nte);

try
    te(1)=readsPar([dname,'.pro'],'alFree[13]');
catch
    te(1)=0;
end

for i=1:nte-1
    te(i+1)=readsPar([dname,'.pro'],sprintf('alFree[%d]',i+13));
end

npez=readsPar([dname,'.pro'],'alFree[11]');

%%
if isold
    
    d=reshape(d,[sz(1:2),ns,npez,sz(3)/npez/ns]);
    d=interleave2linear(d,3);
    d=squeeze(mean(d,4));
    
    
    save([dname,'_zshim'],'d','te');
    
else
    d2=[];
    d1=[];
    
    for i=1:length(com)
        if isempty(com{i})
            d1(:,:,end+1)=d(:,:,i);
        else
            d2(:,:,end+1)=d(:,:,i);
        end
    end
    
    
    d1=d1(:,:,2:end);
    d2=d2(:,:,2:end);
    
    d1=reshape(d1,[sz(1:2),ns,npez,sz(3)/ns/2/npez]);
    d1=squeeze(d1(:,:,:,1,:));
    
    d2=reshape(d2,[sz(1:2),ns,npez,sz(3)/ns/2/npez]);
    d2=squeeze(mean(d2,4));
     
    %d1=interleave2linear(d1);
    %d2=interleave2linear(d2);
    
  %  T2_nozshim=T2map_LinearFit_data(d1,te);
    
  %  T2_zshim=T2map_LinearFit_data(d2,te);
    save([dname,'_NoZshim'],'d1','te');
     
    save([dname,'_zshim'],'d2','te');
    
end


