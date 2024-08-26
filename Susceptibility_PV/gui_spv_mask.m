function p=gui_spv_mask(nTE)

if ~exist('nTE','var')
    nTE = 2;
end

p = parameter('Path2Mask GUI');

p=add(p,'int','swi dir number',12);

%p = add(p,'directoryname',sprintf('swi dcm dir'),'SWI_IMAGES_00??'); 
%p = add(p,'directoryname',sprintf('mag dcm dir'),'MAG_IMAGES_00??'); 
%p = add(p,'directoryname',sprintf('ph dcm dir'),'PHA_IMAGES_00??'); 

p = add(p,'button',sprintf('Sep to %d TE',nTE),'gui_spv_callback(params)'); % separate into two dir; one for each TE

p = add(p,'button','mkdir swi_00#_ana','gui_spv_callback(params)');

p = add(p,'button','Note: Run above in subject folder; below in swi_00#_ana',''); % separate into two dir; one for each TE

p = add(p,'button','Note: load TE2 swi data to image J (file/import/image sequence...)','');

p = add(p,'button','Note: save as analyze in  swi_00#_ana folder','');

p=add(p,'filename','swi analyze file','SWI_IMAGES_00#_TE2.hdr');
p = add(p,'button','invert image','gui_spv_callback(params)'); % separate into two dir; one for each TE
p = add(p,'button','Note: Use inverted images to draw paths in imageJ with neuriteTracer',''); % separate into two dir; one for each TE

%p = add(p,'directoryname','dcm swi dir','../SWI_IMAGES_00#/TE2'); 
p = add(p,'directoryname','neurite path dir','neuriteTracer'); 
p = add(p,'button','path2mask enlarge','gui_spv_callback(params)');
p = add(p,'button','path2mask','gui_spv_callback(params)');

p = add(p,'button','Close','close');

if nargout==0
parametergui(p);
end