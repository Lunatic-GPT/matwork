function knew = kdata_trans_rot(kpos,kdata,shft,rot)
% calculate the k-space data for an object with shift (shft), then rotation
% (rot) given the k-space data before the position change.
% kpos: in units of radian/mm; n*3
% kdata: n*1 or n*1
% shft: in units of mm; 1*3
% rot: in units of degree; [theta, phi, alpha], where alpha is the angle of
% rotation around the axis specified by theta and phi.


shft=repmat(shft,[n,1]);
knew = kdata.*exp(1i*sum(kpos.*shft,2));

rot = 






