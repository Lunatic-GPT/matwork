sid=[3:6,8:14,16:20,22,23,25];
sidc=[4:8,10:23];


for i=1:19
    for j=1:2
        if j==1
            s=sprintf('DB%02d',sid(i));
        else
            s=sprintf('DBControl%02d',sidc(i));
        end



        for k=1:3
            if strcmp(s,'DBControl10') && k>1
                wm=sprintf('../WMmask_shifted/WMmask_%s_scan2.nii.gz',s);
            else
                wm=sprintf('../WMmask_shifted/WMmask_%s.nii.gz',s);
            end

            paname=sprintf('mask/mask_Model_M2EDN2d_nfilter64_61sub_fpm_newwm_Mag_%s_%d.nii.gz',s,k);
            pa=ri(paname);

            wm=ri(wm);

            pa2=pa&wm>0;
            out=sprintf('PAmask_in_wm_%s_%d',s,k);
            save_untouch_niigz(paname,out,pa2);

        end

    end
end
