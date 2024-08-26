function d=detectPVS2SigAllPoints(pc,sig,neg_phase,m_wm)
% m2=detectPVS2SigAllPoints(pc,sig,neg_phase,m_wm)
% detect pixels 2 sigma above zero in zero at all time points.
%


m_wm=ri(m_wm);

a=ri(pc);

if exist('neg_phase','var') && ~isempty(neg_phase) && neg_phase
    a=-a;
else
    neg_phase=false;
end


d=~any(a<=sig*2,4)&m_wm;

tmp=clusterize2(d);
fprintf('%d clusters detected\n',max(tmp(:)));

[dir_name,fname]=fileparts(pc);

prefix=strtok(fname,'.');

save_name=sprintf('%s_VSg%3.1fAll.mat',prefix,sig);
try
    voxsize=ri(m_file,[],[],'vxosize');
    
    center=ri(m_file,[],[],'center');
    
    save(save_name,'d','voxsize','center');
    
catch
    save(save_name,'d');
      
end
