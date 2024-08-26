function tsroi_avetrl(flist,mask,mname,nTR_trial,exc_list)
% tsroi_avetrl(flist,mask,mname,nTR_trial,exc_list)
% mean time course averaged over roi and trials
% flist contain a list of sorted time courses. one for each stimulus pattern.
% exc_list: an cell array of trial indices to be excluded (1 based)
npat = length(flist);
temp = textscan(mask,'%s','Delimiter',';');
mask = maskcalc(temp{1}{:},1,mname);

for i=1:npat
    
[n,d] = nv_roi(flist{i},mask);

if n>0


  d = reshape(d,nTR_trial,length(d)/nTR_trial);
  if exist('exc_list','var')
    d(:,exc_list{i}) = [];
  end
  mn = mean(d,2);
  ts_mn(:,i) = mn - mean(mn(end-4:end));
  ts_sem(:,i) = std(d,0,2)/sqrt(size(d,2));
  ts{i} = d;
end


end

if n>0
    novox = 0;
    save(sprintf('ts_%s.mat',mname),'ts_mn','ts_sem','ts','mask','exc_list','flist','nTR_trial','novox','n');
else
    novox = 1;
    save(sprintf('ts_%s.mat',mname),'novox');
end
