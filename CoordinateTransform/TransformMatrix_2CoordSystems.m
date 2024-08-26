function T=TransformMatrix_2CoordSystems(xyz,xyznew)
% Given the x,y,z axes as the row vecters in xyz ans xyznew, 
% calculate the T matrix that relates v_new=T*v, for vectors in the two
% coord systems


T=xyznew'*xyz;




