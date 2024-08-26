function [nch,chscale]=nChan(scand)

acqmod=readbPar(fullfile(scand,'acqp'),'ACQ_experiment_mode',false);

if ~strcmp('SingleExperiment',acqmod)

    rcvrs=readbPar(fullfile(scand,'acqp'),'ACQ_ReceiverSelect',false);
    nyes=strmatch('Yes',rcvrs);
    nch=length(nyes);    
    chscale=[1.0,0.266,0.1879,0.1769];
else
    nch=1;
    chscale=1;
end
