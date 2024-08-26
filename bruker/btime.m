function t0=btime(id)
%t0=btime(id)
% time in seconds

t0=readbPar(fullfile(id,'acqp'),'ACQ_abs_time');