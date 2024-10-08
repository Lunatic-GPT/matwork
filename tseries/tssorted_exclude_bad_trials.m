function tssorted_exclude_bad_trials(fname_sorted,ts1d,nTR_trl,pk,bl,pk_thr,std_thr)
% tssorted_
% tssorted_exclude_bad_trials(fname_sorted,ts1d,pk,bl,nTR_trl,pk_thr,std_thr)
% 

[d,info] = BrikLoad(fname_sorted);
ntrl = size(d,4)/nTR_trl;

disp('Trials removed:');
ind = [];

if isa(ts1d,'double')
    
  for i=1:ntrl
      if ~any(i==ts1d)
         ind = [ind,(i-1)*nTR_trl+1:i*nTR_trl];
      else
        fprintf('%d, ',i);  
      end
  end
  
  fprintf('\n');

  history = sprintf('tssorted_exclude_bad_trials(%s,ts1d = %s,nTR_trl = %d)',fname_sorted,num2str(ts1d),nTR_trl);

else
    
  for i=1:ntrl
     
    if ~remove_trial(ts1d((i-1)*nTR_trl+1:i*nTR_trl),pk,bl,pk_thr,std_thr)
      ind = [ind,(i-1)*nTR_trl+1:i*nTR_trl];
    else
      fprintf('%d, ',i);
    end
  end

  fprintf('\n');

  history = sprintf('tssorted_exclude_bad_trials(%s,ts1d,pk=%s,bl = %s,nTR_trl = %d,pk_thr = %3.2f,std_thr = %3.2f)',fname_sorted,num2str(pk),num2str(bl),nTR_trl,pk_thr,std_thr);

end

d2 = d(:,:,:,ind);

prefix = strtok(fname_sorted,'+');

WriteBrikEZ(d2,info,history,['rn',prefix,]);

function rm = remove_trial(ts,pk,bl,pk_thr,std_thr)

rm = false;

if std(ts(bl))>std_thr || mean(ts(pk)) - mean(ts(bl)) < pk_thr
    rm = true;
end


