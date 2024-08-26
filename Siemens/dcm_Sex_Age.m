function [sex,age]=dcm_Sex_Age(scand)

age=readdPar(scand,'PatientAge');

sex=readdPar(scand,'PatientSex');

age=str2num(age(1:end-1));
disp(age);