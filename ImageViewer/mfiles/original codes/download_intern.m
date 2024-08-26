function download_intern

load download_image_gui 

DCE_no = get(params,'DCE #');
DTI_no = get(params,'DTI #');

fnumber_array = [DCE_no DTI_no];
%sourcy = get(params, 'Enter the file name');
%sourcy = [sourcy '/'];
sourcy2 = get(params, 'Select source directory name');
dest2 = get(params, 'Select destination directory name');
DTI_iter = get(params, 'DTI iteration');
DCE_iter = get(params, 'DCE iteration');

% delta1 = get(params, 'delta1 [sec]');
% delta2 = get(params, 'delta2 [sec]');
% delta3 = get(params, 'delta3 [sec]');
% delta4 = get(params, 'delta4 [sec]');
% delta5 = get(params, 'delta5 [sec]');
% delta_vec = [delta1 delta2 delta3 delta4 delta5];
for i=1:4
    if i==1
        type = 'DCE';
    elseif i ==2
        type = 'ADC';
    elseif i == 3
        type = 'DTI';
    else
        type = 'FA';
    end
    source = [sourcy2 '/' type];
    dest = [dest2 '/'];
    if i == 1
        fnumber = fnumber_array(1);
    else
        fnumber = fnumber_array(2);
    end
    recon_hahn(i,1, source,DTI_iter,dest,fnumber, DCE_iter);
end

%close all;