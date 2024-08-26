function set_tpwrtag(nTR_pre,nTR_trl,pwr,fname)
% set_tpwrtag(nTR_pre,nTR_trl,pwr,fname)
% set the pwr levels to pwr array values 
% pwr: an array of power levels for the trials.
t = zeros(nTR_pre,1);
for i=1:length(pwr)
 t(end+1:end+nTR_trl) = pwr(i);
end

t(end+1) = 0;

save_mat_int(t,fname);