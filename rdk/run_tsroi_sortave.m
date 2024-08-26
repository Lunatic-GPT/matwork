function run_tsroi_sortave(p_ts,mstr,mname)
% replace old run_tsroi_arearoi_sortave 5-7-2010
if ~exist('results','dir')
    mkdir('results');
end

p_ts = set(p_ts,'roi mask',mstr);
fname = sprintf('results/ts_%s',mname);
p_ts=set(p_ts,'save file name',[fname,'.mat']);
tsroi_sortave(p_ts);
saveas(gcf,[fname,'.fig']);