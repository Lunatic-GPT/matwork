function CBVa_MT(flist_dtr,flist,bl,prefix_out)
% CBVa_MT(flist_dtr,flist,bl,prefix_out)
% the files should be ordered in the descending order of tpwrtag.

for i=1:length(flist)
    if i==1
        [d,info]=BrikLoad(flist{i});
        d_dtr = BrikLoad(flist_dtr{i});
    else
        tmp = BrikLoad(flist{i});
        d=cat(5,d,tmp);
        tmp = BrikLoad(flist_dtr{i});
        d_dtr=cat(5,d_dtr,tmp);
    end
end

S0=mean(d(:,:,:,bl,:),4);


MTR = 1-S0./repmat(S0(:,:,:,:,end),[1,1,1,1,size(d,5)]);
MTR = permute(MTR,[1,2,3,5,4]);
sz = size(d_dtr);
ds = d_dtr-repmat(mean(d_dtr(:,:,:,bl,:),4),[1,1,1,sz(4),1]);
ds_norm = ds./repmat(S0(:,:,:,:,end),[1,1,1,sz(4),sz(5)]);
ds_norm_BOLD = ds./repmat(S0,[1,1,1,sz(4),1]);  
tmp = S0(:,:,:,1,end);
tmp = sort(tmp(:));
imax = mean(tmp(end-100:end));

bd = zeros(sz(1:4));
tissue = zeros(sz(1:4));  %relative signal change in tissue.
for i=1:sz(1)
    for j=1:sz(2)
        for k=1:sz(3)
           % if S0(i,j,k,1,end)<imax/10
            %    continue;
            %end
            
            for l=1:sz(4)
             x=1-MTR(i,j,k,:);
             y = squeeze(ds_norm(i,j,k,l,:));
             b=lscov([x(:),ones(size(x(:)))],y(:));
             bd(i,j,k,l) = b(2);
             tissue(i,j,k,l) = b(1);
            end
        end
    end
    disp(i);
end

str2=sprintf('%s;',flist{:});

str1=sprintf('%s;',flist_dtr{:});
history = sprintf('CBVa_MT(%s,%s,%s,%s)',str1,str2,num2str(bl),prefix_out);
WriteBrikEZ(bd,info,history,prefix_out);
WriteBrikEZ(tissue,info,history,[prefix_out,'_slope']);

WriteBrikEZ(MTR,info,history,[prefix_out,'_MTR']);

for i=1:length(flist)
    prefix = strtok(flist_dtr{i},'+');
    WriteBrikEZ(ds_norm(:,:,:,:,i),info,history,[prefix,'_normMT0']);
    WriteBrikEZ(ds_norm_BOLD(:,:,:,:,i),info,history,[prefix,'_norm']);
end




