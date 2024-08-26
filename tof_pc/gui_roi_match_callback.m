function gui_roi_match_callback(param)

m1=get_fpattern(param,'mask 1');
m2=get_fpattern(param,'mask 2');

u1=get_fpattern(param,'underlay 1');
u2=get_fpattern(param,'underlay 2');

max_dist=get(param,'max dist (mm)');
vox_size=get(param,'vox size (mm)');

d1=ri_d1(m1);
d2=ri_d1(m2);

if ~isempty(u1)
 u1=ri(u1);
 u2=ri(u2);
end

[d1n,d2n,ind1,ind2]=clusters_match(d1,d2,round(max_dist/vox_size),false);  

save(filename_append(m1,'_match_ind'),'ind1','ind2');


disp('gui_cm_callback done!');



