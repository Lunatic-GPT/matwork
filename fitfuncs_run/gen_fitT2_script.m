%fid=fopen('run_fit_T2_ase_nozshim.sh','w');
fid=fopen('run_fit_T2_slicebyslice.sh','w');

for i=1:20
    
fprintf(fid,'bsub matbgk \"fit_T2_ase_zshim(%d:%d)\" log%d\n',i*2-1,i*2,i);

end

fclose(fid);

