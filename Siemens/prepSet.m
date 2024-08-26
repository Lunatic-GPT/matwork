function [m_dPhase,m_dFrequency]=prepSet(gro,tro,SliceOffCenter,FOV,thickness,lpe,l3d)

%gro in mT/m
%tro: the readout time (us);
%sliceoffcenter: 1*3 (mm)(pe,ro,3d)% +: (in the directon defined by the first, 2nd, and 3rd columns of Rotmat?) 

%where Rotmat: x: R->L; y: P->A; or z: S->I;

% FOV: 1*3 (mm) (pe,ro,3d); include oversampling along pe and 3d
% thickness: slice thickness (mm);
%lpe, l3d: pe and 3D indices, center 0;

o3d=SliceOffCenter(3); %
ope=SliceOffCenter(1); %pe
oro=SliceOffCenter(2); %ro


PhaseOffCenter3D = 360*(thickness/2+o3d)/FOV(3);  % parameters in sSLICE_POS
PhaseOffCenterPE = 360*(ope)/FOV(1);

disp([PhaseOffCenterPE,PhaseOffCenter3D]);
m_dFrequency = gro*oro*0.001*42580;%Hz
m_dPhase = -360*tro/2*m_dFrequency/1000000+PhaseOffCenterPE*lpe+PhaseOffCenter3D*l3d;
m_dPhaseNeg = -(m_dPhase+m_dFrequency*tro/1000000*360);

disp([m_dFrequency,m_dPhase]);
disp(m_dPhaseNeg);

