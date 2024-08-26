function p=gui_veinMoment

p = parameter('Moment Vein GUI');

p = add(p,'directoryname',sprintf('mag dcm dir'),'../MAG_IMAGES_00??/TE2'); 
p = add(p,'directoryname',sprintf('ph dcm dir'),'../PHA_IMAGES_00??/TE2'); 
p=add(p,'float','B0 (T)',7);

p = add(p,'filename','wm mask for Frangi','mask_wm.nii'); 

p = add(p,'button','calc moment (Frangi)','gui_veinMoment_callback(params)'); % separate into two dir; one for each TE


%p = add(p,'directoryname','neurite path dir','neuriteTracr'); 
%p = add(p,'string','prefix','neuriteTracer'); 

%p = add(p,'button','path to mask','gui_veinMoment_callback(params)');
p=add(p,'filename','path2mask output','mask_neuriteTrace.mat');
p = add(p,'button','calc_moment_Neurite','gui_veinMoment_callback(params)'); %
p = add(p,'button','gen bsub script','gui_veinMoment_callback(params)'); 
p=add(p,'float','dchi (ppm)',0.45);
p = add(p,'button','show histogram','gui_veinMoment_callback(params)'); 

p = add(p,'button','Close','close');

if nargout==0
parametergui(p);
end


