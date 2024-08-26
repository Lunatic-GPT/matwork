function p=gui_resliceT1()

p = parameter('reslice T1');

%p=add(p,'int','Files',nraw);
%for i=1:nraw
p=add(p,'directoryname','T1 image','GRE_IR');  % leave blank if use the fid file in scan folder.

p = add(p,'int','center shift','0 0 0');

p = add(p,'float','new dimension (mm)','0.1563 0.1563 2');
p = add(p,'float','new matrix','1280 1040 1');

p = add(p,'int','flip xy',1);

p = add(p,'int','flip ap',1);

p = add(p,'int','flip lr',0);

p=add(p,'button','reslice','gui_resliceT1_callback(params)');

p = add(p,'button','Close','close');

if nargout==0
p = parametergui(p);
end




