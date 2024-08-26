function save_mc_data(h5file, mcfile)
%mcfile should contain data named reskESPIRiT
a=h5read([prefix,'.h5'],'/dataset/data');
load(mcfile);
%%
data_save=[]; %nlines*(nro*nch)

%zeros(nro,nch)
for i=1:length(a.data)
    
    i1=a.head.idx.kspace_encode_step_1(i)+1;
    i2=a.head.idx.kspace_encode_step_2(i)+1;
    flags=a.head.flags(i); 

    if flags~=8388608 && flags~=262144 &&flags~=134217728 %noise scan; first two repetitions at k=0; fatnav 

       tmp=squeeze(reskESPIRiT(:,i1,i2,:)); 
       data_save=cat(1,data_save,transpose(tmp(:)));
        fprintf('%d %d %d  %d  %d\n',i,i1,i2,flags,length(a.data{i}));
    end

end
%%
[d,~,~,feed_back]=Read_UIH_Raw_v5_7_additional_shot_index('UID_7253832852478711979_fse_mx_3d_fatnav2_ACS5.raw');

d=squeeze(d);
d2=sos(sos(d,3),4);

%d: 352*270*704*64 (ro in the 3rd dimension):
%terms in pe par plane;
% n_32768 = 1907; the difference in the k-space center

%feed_back: 148*64*68*58*58; 
%data from two extra repetitions: total 8 repetitions 


