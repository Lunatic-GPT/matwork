function res=FOVcenter_dcm2pro(a)

if ~exist([a,'.pro'],'file')   
 extp(a);
end

 posx=readsPar([a,'.pro'],'asSlice[0].sPosition.dSag');
  posy=readsPar([a,'.pro'],'asSlice[0].sPosition.dCor');
  posz=readsPar([a,'.pro'],'asSlice[0].sPosition.dTra');
     
res=[posx,posy,posz];