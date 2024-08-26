function renbScans(maxScan)

if ~exist('maxScan','var')
    maxScan=150;
end

fid=fopen('renbScans.log','w');

for i=1:maxScan
    
    orig=num2str(i);
    if exist(num2str(i),'dir')
     pulprog=readbPar(fullfile(num2str(i),'acqp'),'ACQ_protocol_name',false);
     pulprog=pulprog{1};
    end
    
    if exist(fullfile(orig,'fid'),'file')
     nf=dir([pulprog(2:end-1),'*']);
     dest=sprintf('%s_%d',pulprog(2:end-1),length(nf)+1);
     if ~exist(dest,'dir')
       fprintf('Rename %s to %s\n',orig,dest);
       fprintf(fid,'Rename %s to %s\n',orig,dest);
       
       movefile(orig,dest,'f');
     end
    end
end

fclose(fid);