% s = zeros(512,512);
test = cat(3,image3, image4, image5, image6);
for i = 1: 512
    for j = 1:512
        temp = (squeeze(test(i,j,:)))';
        fit_vector=polyfit([1 2 3 4],temp,1);
        slope = fit_vector(1);
        s(i,j) = slope;
    end
end
save('C:\Documents and Settings\hahntobi\My Documents\MATLAB\temporary saved variables\wash_out', 's');