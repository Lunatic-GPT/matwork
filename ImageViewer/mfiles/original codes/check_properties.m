function check_properties
% check_properties loads images from the given data files and calculates
% the given parameters

sourcy = input('Please enter the patient ID: ','s'); %% CHANGE THIS
sourcy1 = [sourcy '/'];
sourcy2 = ['~/matlab/patient_images/' sourcy1 '/DCE_temp/'];
im_mean = zeros(18,1);
for i = 1:9
    number = num2str(i);
    filename = ['0' number];
    source = [sourcy2 filename '.mat'];
    load(source, 'image');
    im_mean(2*i-2) = mean(mean(image));
    im_mean(2*i-1) = max(max(image));
end
for i=10:18
    filename = num2str(i);
    source = [sourcy2 filename '.mat'];
    load(source, 'image');
    im_mean(i) = mean(mean(image));
end

disp('The mean vector min and max')
disp(floor(min(im_mean)));
disp(ceil(max(im_mean)));


    