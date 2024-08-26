!@auto_tlrc -apar anat+tlrc -input V5u+orig -dxyz 3 -rmode NN
!@auto_tlrc -apar anat+tlrc -input V5m+orig -dxyz 3 -rmode NN
for i=1:3
    cmd = sprintf('@auto_tlrc -apar anat+tlrc -input ts_sort_coh%d+orig -dxyz 3',i);
    unix(cmd);
end
