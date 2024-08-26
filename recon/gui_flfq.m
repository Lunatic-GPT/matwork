function gui_flfq()

p = parameter('gui to recon fl fq data');

%p=add(p,'int','Files',nraw);
%for i=1:nraw
p=add(p,'filename',sprintf('raw data'),'*_sb_*.dat');  % leave blank if use the fid file in scan folder.

p = add(p,'button','raw to mat','gui_flfq_callback(params)');
%p=add(p,'button','FatNav recon script','gui_flfq_callback(params)');

%end
p=add(p,'int','cardiac phases',12);

%% this two files should be located in the folders named as the raw data prefix.
% ok if do not exist.

p=add(p,'filename','real peaks',''); 
p=add(p,'filename','fake peaks','');
p=add(p,'filename','MON peaks','');  % peaks that should be replaced by the middle point of two neighboring peaks
p=add(p,'filename','interp n peaks','');  %add two peaks between a pair of peaks a third number specifies the number of peaks to add.


p=add(p,'float','max heart rate (per min)',80);

p = add(p,'button','load physio','gui_flfq_callback(params)');
p = add(p,'button','check pulse/gen bsub','gui_flfq_callback(params)');


p=add(p,'float','interp',[2,3.3333]);  % factor after keeping only the nonzero fraction

%p=add(p,'float','nonzero fraction',1);

p=add(p,'button','NOTE: to run bsub, upload raw or mat data file and its folder containing','');
p=add(p,'button','.pro files to the server, also the respiratory peak txt files.','');
p=add(p,'button','gen recon bsub','gui_flfq_callback(params)');
p=add(p,'filename','SB ref','');

p=add(p,'button','recon','gui_flfq_callback(params)');

p=add(p,'button','gen mb recon bsub','gui_flfq_callback(params)');
p = add(p,'button','Clear','');

p = add(p,'button','Close','close');

p = parametergui(p);



