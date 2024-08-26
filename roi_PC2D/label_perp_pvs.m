function m2=label_perp_pvs(fname,ind_pvs,prefix)
% label ind_pvs to 2, others to 1;

m=ri_d1(fname);

m2=int32(m>0);
for i=1:length(ind_pvs)
   
    m2(ind_pvs(i)==m)=2;
    
end

save(prefix,'m2');

