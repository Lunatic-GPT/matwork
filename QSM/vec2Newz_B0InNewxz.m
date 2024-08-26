function res = vec2Newz_B0InNewxz(vec,B0)

[th,phi]=unitVec2thetaPhi(vec(:));

m1=transform_matrix_rotation(0,-phi);

m2=transform_matrix_rotation_arb_axis(90,90,-th);

res=m2*m1;

B0p=res*B0(:);

[th2,phi2]=unitVec2thetaPhi(B0p(:));

m3=transform_matrix_rotation(0,-phi2);

res=m3*res;

