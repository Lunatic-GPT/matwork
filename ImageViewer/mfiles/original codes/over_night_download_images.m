function over_night_download_images
% over_night_download_images constructs all the given images, calculates
% and saves the wash-in and wash-out images etc.

close all;
disp('1) DCE, 2) ADC, 3) DTI, 4) FA: ');

%% FNUMBER_ARRAY: [DCE DTI]
% 
% fnumber_array = [106 21]; %% CHANGE THIS
% sourcy = 'N84432/'; %% CHANGE THIS
% DTI_iter = 33; % CHANGE THIS
% DCE_iter = 224;
% % delta_vec = [223 92 91 92 92];
% for i=1:4
%     if i==1
%         type = 'DCE';
%     elseif i ==2
%         type = 'ADC';
%     elseif i == 3
%         type = 'DTI';
%     else
%         type = 'FA';
%     end
%     source = ['/export/res/breast/' sourcy type];
%     dest = ['~/matlab/patient_images/' sourcy];
%     if i == 1
%         fnumber = fnumber_array(1);
%     else
%         fnumber = fnumber_array(2);
%     end
%     recon_hahn(i,2, source,DTI_iter,dest,fnumber, DCE_iter,8);
% end
% 
% fnumber_array = [125 20]; %% CHANGE THIS
% sourcy = 'M179579/'; %% CHANGE THIS
% DTI_iter = 33; % CHANGE THIS
% DCE_iter = 224;
% % delta_vec = [203 258 92 91 92];
% for i=1:4
%     if i==1
%         type = 'DCE';
%     elseif i ==2
%         type = 'ADC';
%     elseif i == 3
%         type = 'DTI';
%     else
%         type = 'FA';
%     end
%     source = ['/export/res/breast/' sourcy type];
%     dest = ['~/matlab/patient_images/' sourcy];
%     if i == 1
%         fnumber = fnumber_array(1);
%     else
%         fnumber = fnumber_array(2);
%     end
%     recon_hahn(i,2, source,DTI_iter,dest,fnumber, DCE_iter,12);
% end
% 
% 
fnumber_array = [101 18];
sourcy = 'M298821/'; %% CHANGE THIS
DTI_iter = 33; % CHANGE THIS
DCE_iter = 224;
% delta_vec = [361 223 90 91 91];
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
    source = ['/export/res/breast/' sourcy type];
    dest = ['~/matlab/patient_images/' sourcy];
    if i == 1
        fnumber = fnumber_array(1);
    else
        fnumber = fnumber_array(2);
    end
    recon_hahn(i,2, source,DTI_iter,dest,fnumber, DCE_iter,29);
end

% fnumber_array = [125 21]; %% CHANGE THIS
% sourcy = 'M297672/'; %% CHANGE THIS
% DTI_iter = 33; % CHANGE THIS
% DCE_iter = 224;
% % delta_vec = [278 196 90 91 91];
% for i=1:4
%     if i==1
%         type = 'DCE';
%     elseif i ==2
%         type = 'ADC';
%     elseif i == 3
%         type = 'DTI';
%     else
%         type = 'FA';
%     end
%     source = ['/export/res/breast/' sourcy type];
%     dest = ['~/matlab/patient_images/' sourcy];
%     if i == 1
%         fnumber = fnumber_array(1);
%     else
%         fnumber = fnumber_array(2);
%     end
%     recon_hahn(i,2, source,DTI_iter,dest,fnumber, DCE_iter,12);
% end
% 
% 
% fnumber_array = [112 20]; %% CHANGE THIS
% sourcy = 'B299986/'; %% CHANGE THIS
% DTI_iter = 33; % CHANGE THIS
% DCE_iter = 224;
% % delta_vec = [133 189 91 91 91];
% for i=1:4
%     if i==1
%         type = 'DCE';
%     elseif i ==2
%         type = 'ADC';
%     elseif i == 3
%         type = 'DTI';
%     else
%         type = 'FA';
%     end
%     source = ['/export/res/breast/' sourcy type];
%     dest = ['~/matlab/patient_images/' sourcy];
%     if i == 1
%         fnumber = fnumber_array(1);
%     else
%         fnumber = fnumber_array(2);
%     end
%     recon_hahn(i,2, source,DTI_iter,dest,fnumber, DCE_iter,10);
% end
% 
% fnumber_array = [121 18]; %% CHANGE THIS
% sourcy = 'B298964/'; %% CHANGE THIS
% DTI_iter = 33; % CHANGE THIS
% DCE_iter = 224;
% % delta_vec = [317 240 91 91 91];
% for i=1:4
%     if i==1
%         type = 'DCE';
%     elseif i ==2
%         type = 'ADC';
%     elseif i == 3
%         type = 'DTI';
%     else
%         type = 'FA';
%     end
%     source = ['/export/res/breast/' sourcy type];
%     dest = ['~/matlab/patient_images/' sourcy];
%     if i == 1
%         fnumber = fnumber_array(1);
%     else
%         fnumber = fnumber_array(2);
%     end
%     recon_hahn(i,2, source,DTI_iter,dest,fnumber, DCE_iter, 9);
% end
% 
% fnumber_array = [83 12]; %% CHANGE THIS
% sourcy = 'B153241/'; %% CHANGE THIS
% DTI_iter = 33; % CHANGE THIS
% DCE_iter = 224;
% % delta_vec = [307 181 91 91 91];
% for i=1:4
%     if i==1
%         type = 'DCE';
%     elseif i ==2
%         type = 'ADC';
%     elseif i == 3
%         type = 'DTI';
%     else
%         type = 'FA';
%     end
%     source = ['/export/res/breast/' sourcy type];
%     dest = ['~/matlab/patient_images/' sourcy];
%     if i == 1
%         fnumber = fnumber_array(1);
%     else
%         fnumber = fnumber_array(2);
%     end
%     recon_hahn(i,2, source,DTI_iter,dest,fnumber, DCE_iter,5);
% end

fnumber_array = [106 18]; %% CHANGE THIS
sourcy = 'B162582/'; %% CHANGE THIS
DTI_iter = 33; % CHANGE THIS
DCE_iter = 224;
% delta_vec = [263 200 92 92 92];
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
    source = ['/export/res/breast/' sourcy type];
    dest = ['~/matlab/patient_images/' sourcy];
    if i == 1
        fnumber = fnumber_array(1);
    else
        fnumber = fnumber_array(2);
    end
    recon_hahn(i,2, source,DTI_iter,dest,fnumber, DCE_iter,17);
end

% fnumber_array = [65 6]; %% CHANGE THIS
% sourcy = 'M273920/'; %% CHANGE THIS
% DTI_iter = 33; % CHANGE THIS
% DCE_iter = 224;
% % delta_vec = [317 219 90 90 90];
% for i=1:4
%     if i==1
%         type = 'DCE';
%     elseif i ==2
%         type = 'ADC';
%     elseif i == 3
%         type = 'DTI';
%     else
%         type = 'FA';
%     end
%     source = ['/export/res/breast/' sourcy type];
%     dest = ['~/matlab/patient_images/' sourcy];
%     if i == 1
%         fnumber = fnumber_array(1);
%     else
%         fnumber = fnumber_array(2);
%     end
%     recon_hahn(i,2, source,DTI_iter,dest,fnumber, DCE_iter,8);
% end
% 
% fnumber_array = [88 15]; %% CHANGE THIS
% sourcy = 'M301431/'; %% CHANGE THIS
% DTI_iter = 33; % CHANGE THIS
% DCE_iter = 224;
% % delta_vec = [297 269 92 91 92];
% for i=1:4
%     if i==1
%         type = 'DCE';
%     elseif i ==2
%         type = 'ADC';
%     elseif i == 3
%         type = 'DTI';
%     else
%         type = 'FA';
%     end
%     source = ['/export/res/breast/' sourcy type];
%     dest = ['~/matlab/patient_images/' sourcy];
%     if i == 1
%         fnumber = fnumber_array(1);
%     else
%         fnumber = fnumber_array(2);
%     end
%     recon_hahn(i,2, source,DTI_iter,dest,fnumber, DCE_iter,19);
% end
% % 
% 
% fnumber_array = [63 11]; %% CHANGE THIS
% sourcy = 'B49217/'; %% CHANGE THIS
% DTI_iter = 35; % CHANGE THIS
% DCE_iter = 224;
% % delta_vec = [200 243 92 92 91];
% for i=1:4
%     if i==1
%         type = 'DCE';
%     elseif i ==2
%         type = 'ADC';
%     elseif i == 3
%         type = 'DTI';
%     else
%         type = 'FA';
%     end
%     source = ['/export/res/breast/' sourcy type];
%     dest = ['~/matlab/patient_images/' sourcy];
%     if i == 1
%         fnumber = fnumber_array(1);
%     else
%         fnumber = fnumber_array(2);
%     end
%     recon_hahn(i,2, source,DTI_iter,dest,fnumber, DCE_iter,5);
% end
% % 

fnumber_array = [112 22];
sourcy = 'M303424/'; %% CHANGE THIS
DTI_iter = 33; % CHANGE THIS
DCE_iter = 224;
% delta_vec = [276 194 92 91 92];
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
    source = ['/export/res/breast/' sourcy type];
    dest = ['~/matlab/patient_images/' sourcy];
    if i == 1
        fnumber = fnumber_array(1);
    else
        fnumber = fnumber_array(2);
    end
    recon_hahn(i,2, source,DTI_iter,dest,fnumber, DCE_iter,16);
end
% 
% fnumber_array = [114 11];
% sourcy = 'M7803/'; %% CHANGE THIS
% DTI_iter = 33; % CHANGE THIS
% DCE_iter = 224;
% % delta_vec = [226 268 91 92 92];
% for i=1:4
%     if i==1
%         type = 'DCE';
%     elseif i ==2
%         type = 'ADC';
%     elseif i == 3
%         type = 'DTI';
%     else
%         type = 'FA';
%     end
%     source = ['/export/res/breast/' sourcy type];
%     dest = ['~/matlab/patient_images/' sourcy];
%     if i == 1
%         fnumber = fnumber_array(1);
%     else
%         fnumber = fnumber_array(2);
%     end
%     recon_hahn(i,2, source,DTI_iter,dest,fnumber, DCE_iter,21);
% end
% 
% fnumber_array = [67 11]; %% CHANGE THIS
% sourcy = 'B49217/'; %% CHANGE THIS
% DTI_iter = 35; % CHANGE THIS
% DCE_iter = 224;
% % delta_vec = [200 243 92 92 91];
% for i=1:4
%     if i==1
%         type = 'DCE';
%     elseif i ==2
%         type = 'ADC';
%     elseif i == 3
%         type = 'DTI';
%     else
%         type = 'FA';
%     end
%     source = ['/export/res/breast/' sourcy type];
%     dest = ['~/matlab/patient_images/' 'B49217b/'];
%     if i == 1
%         fnumber = fnumber_array(1);
%     else
%         fnumber = fnumber_array(2);
%     end
%     recon_hahn(i,2, source,DTI_iter,dest,fnumber, DCE_iter,7);
% end
% 
% fnumber_array = [120 11]; 
% sourcy = 'N309953/';  % MALIGN
% DTI_iter = 33; %  !!!!!!! THE DTI NUMBERS HERE ARE NOT CORRECT - WE DON'T USE THIS CASE FOR THE ADC (...) - CALCULATION! (BECAUSE OF NEGROSIS...)
% DCE_iter = 224;
% % delta_vec = [200 243 92 92 91];
% for i=1:4
%     if i==1
%         type = 'DCE';
%     elseif i ==2
%         type = 'ADC';
%     elseif i == 3
%         type = 'DTI';
%     else
%         type = 'FA';
%     end
%     source = ['/export/res/breast/' sourcy type];
%     dest = ['~/matlab/patient_images/' sourcy];
%     if i == 1
%         fnumber = fnumber_array(1);
%     else
%         fnumber = fnumber_array(2);
%     end
%     recon_hahn(i,2, source,DTI_iter,dest,fnumber, DCE_iter,11);
% end
% 
% fnumber_array = [110 20]; 
% sourcy = 'N309953/';  % BENIGN
% DTI_iter = 33; 
% DCE_iter = 224;
% for i=1:4
%     if i==1
%         type = 'DCE';
%     elseif i ==2
%         type = 'ADC';
%     elseif i == 3
%         type = 'DTI';
%     else
%         type = 'FA';
%     end
%     source = ['/export/res/breast/' sourcy type];
%     dest = ['~/matlab/patient_images/' 'N309953benign'];
%     if i == 1
%         fnumber = fnumber_array(1);
%     else
%         fnumber = fnumber_array(2);
%     end
%     recon_hahn(i,2, source,DTI_iter,dest,fnumber, DCE_iter,10);
% end