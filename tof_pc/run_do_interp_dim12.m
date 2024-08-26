
dname={'meas_MID98_fl_fq_retroZ_mb_0_3by0_5mm_FID15385'};

for i=1

    dsb=readMeasDat([dname{i},'.dat'],32,0,true);    
    mid=strtok_no(dname{i},'_',2);
    cur_dir= cd(dname{i});
    a=ri(['recon_',mid,'.mat']);
    
    interp_factor=[3,5];
    do_interp_dim12(a,interp_factor,mid,dsb(:,1:32));

    cd(cur_dir);
end

 