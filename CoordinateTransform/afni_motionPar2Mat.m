function m_all=afni_motionPar2Mat(dfile)
% n  roll  pitch  yaw  dS  dL  dP  rmsold rmsnew
%                 roll  = rotation about the I-S axis }
%                    pitch = rotation about the R-L axis } degrees CCW
%                    yaw   = rotation about the A-P axis }
%                      dS  = displacement in the Superior direction  }
%                      dL  = displacement in the Left direction      } mm
%                      dP  = displacement in the Posterior direction }
% so the motion parameters are the amount of motion from new to base.                      
% dicom convention; right-hand coordinate system
% the x-axis is increasing to the left hand side of the patient. 
% The y-axis is increasing to the posterior side of the patient.
% The z-axis is increasing toward the head of the patient.
%#define ORI_R2L_TYPE  0  /* Right to Left         */
%#define ORI_L2R_TYPE  1  /* Left to Right         */
%#define ORI_P2A_TYPE  2  /* Posterior to Anterior */
%#define ORI_A2P_TYPE  3  /* Anterior to Posterior */
%#define ORI_I2S_TYPE  4  /* Inferior to Superior  */
%#define ORI_S2I_TYPE  5  /* Superior to Inferior  */
% the AFNI convention is that R-L, A-P, and I-S are
%         negative-to-positive. same as DICOM.
% RL 1 dim: AP: 2nd dim; IS: 3rd dim
% the output should be reshaped to ([4,3]);
% the output is matrix (m) is defined such that x=m*x', where x and x' are the original and
% new coordinates of the same point of the object. 
if isa(dfile,'char')
d=load(dfile);
else
    d=dfile;
end

d=d(:,2:7);
m_all=zeros(size(d,1),12);
for i=1:size(d,1)

    m_IS = transform_matrix_rotation_arb_axis(0,0,d(i,1));
    m_RL = transform_matrix_rotation_arb_axis(90,0,d(i,2));
    m_AP = transform_matrix_rotation_arb_axis(90,90,d(i,3));
    
    %m=m_IS*m_RL*m_AP;
    %m=m_IS*m_AP*m_RL;
    %m=m_AP*m_RL*m_IS;
    %m=m_AP*m_IS*m_RL;
    %m=m_RL*m_IS*m_AP;
    m=m_AP*m_RL*m_IS;  % not sure which one is correct
    
   % shift=m*d(i,[5,6,4])';
    m(:,end+1)=d(i,[5,6,4]);
    
m_all(i,:)=vec(m');
end  

