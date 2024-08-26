function sbs(maxScan)

if ~exist('maxScan','var')
    maxScan=150;
end


for i=1:maxScan
    
    if exist(num2str(i),'dir') && exist(fullfile(num2str(i),'fid'),'file')
     pulprog=readbPar(fullfile(num2str(i),'acqp'),'ACQ_protocol_name',false);
     pulprog=pulprog{1};
     
     fprintf('%d %s\n',i,pulprog);
    end
    
end
