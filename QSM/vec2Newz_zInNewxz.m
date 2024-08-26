function res = vec2Newz_zInNewxz(vec)

[th,phi]=unitVec2thetaPhi(vec(:));

m1=transform_matrix_rotation(0,-phi);

m2=transform_matrix_rotation_arb_axis(90,90,-th);

res=m2*m1;
