function init_mat_flirt(theta,phi,deg,shift)


m=transform_matrix_rotation_arb_axis(theta,phi,deg);

m(:,4)=shift';
m(4,:)=[0,0,0,1];
disp(m);
