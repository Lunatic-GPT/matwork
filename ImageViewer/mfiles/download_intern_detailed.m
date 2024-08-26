function download_intern_detailed(params)

%load download_image_gui 

set(gcbo,'Enable','off');
sourcy2 = get(params, 'Select source directory name');
lesion = get(params, 'lesion name');
DCE_iter = get(params, 'total slices/volume');
fslice = get(params, 'first slice');
lslice = get(params,'last slice');


dest2 = fullfile(sourcy2,'..',lesion);
if ~exist(dest2,'dir')
    mkdir(dest2);
end

recon_DCE(sourcy2,fullfile(dest2,'DCE_temp'),fslice, DCE_iter, lslice-fslice+1);

save(fullfile(dest2,'DCE_temp','download_image_gui.mat'),'params');

set(gcbo,'Enable','on');
cd(dest2);


%close all;