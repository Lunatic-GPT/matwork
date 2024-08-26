function mask=mask_veinMoment(sz,rad,voxSize,center)

% sz: 1*2;
% rad: the inner and outer radii


i_rad_mom=round(rad/voxSize);

                    
   m_momtmp = mask_circle(sz,i_rad_mom(1),center,1);
   m_momtmp2 = mask_circle(sz,i_rad_mom(2),center,1);
   mask = m_momtmp==0&m_momtmp2>0;
                        