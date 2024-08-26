function xform_new=TransformMatrix_NewCoordSystem(xform_old,T)
% xform_old is the xform in the old system that relates two vectors
% u=xform*v
% xform_new is the xform in the new system that relates the same vectors in
% the new coord system
% u_new=xform_new*v_new;
% T is the transformation that relates the coordiates in the new and old
% systems v_new=T*v and u_new=T*u;

xform_new=T*xform_old*T';




