function p=gui_FatNav()

p = parameter('gui to recon FatNav data');
p=add(p,'string','sid','');
p=add(p,'filename','raw file','');  % leave blank if use the fid file in scan folder.
p=add(p,'int','images per process','20');  % step great than 10 can exceed memory limit.

p=add(p,'button','Extract Data','gui_FatNav_callback(params)');

p=add(p,'button','gen bsub','gui_FatNav_callback(params)');


p=add(p,'button','push to longleaf','gui_FatNav_callback(params)');

p=add(p,'button','run','gui_FatNav_callback(params)');

p=add(p,'button','gen bsub+push to longleaf+run','gui_FatNav_callback(params)');
p=add(p,'button','Extract Data+gen bsub+push to killdevil/longleaf+run','gui_FatNav_callback(params)');

%p=add(p,'button','get output files','gui_FatNav_callback(params)');

%p=add(p,'button','Extract Data+gen bsub+push to killdevil/longleaf+run+get output','gui_FatNav_callback(params)');

p=add(p,'button','Note: run above under raw data folder',''); 
%p=add(p,'int','Repetitions: start stop step','1 128 10');  % step great than 10 can exceed memory limit.
%p=add(p,'button','gen bsub script','gui_FatNav_callback(params)');

p=add(p,'int','Repititions: start stop','1 128');
p=add(p,'button','recon','gui_FatNav_callback(params)');

p=add(p,'int','Repititions: start','');
p=add(p,'button','mat2afni','gui_FatNav_callback(params)');

p=add(p,'button','volreg','gui_FatNav_callback(params)');

 p=add(p,'filename','dfile','Motion_MID.1D');  % leave blank if use the fid file in scan folder.

p=add(p,'button','plot motion','gui_FatNav_callback(params)');
p=add(p,'int','number of map calc jobs',10);
p=add(p,'int','number of nufft blocks',1);  % use 16 for saving memory
p=add(p,'button','recon T2 w/mc bsub','gui_FatNav_callback(params)');

p=add(p,'int','iro indices (start stop step)','1 512 20');
p=add(p,'int','maps',2);

p=add(p,'button','ESPIRiT bsub','gui_FatNav_callback(params)');


%p=add(p,'int','Files',nraw);
%for i=1:nraw

%end

p = add(p,'button','Close','close');
if nargout==0
parametergui(p);
end


