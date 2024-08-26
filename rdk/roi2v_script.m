sid = 'Hua';

roil = {'V1','V2d','V2v','V3d','V3v','V4','MTLoc'};
roir = roil;
suma = '../../segmentation/SUMA';
if strcmp(sid,'BQ')
    roil = {'V1','V2dx','V2v','V3dx','V3v','V4','MTLoc'}; %Peng
    suma = '../../SUMA'; %peng
end
if strcmp(sid,'Hua')
   roil = {'V1','V2d','V2v','V3dx','V3v','V4','MTLoc'}; %Hua
end
if strcmp(sid,'K004')
   roir = {'V1','V2d','V2v','V3d','V3v','V4'};
   roil= roir;
end

%
for i=1:length(roil)
    roil{i} = [roil{i},'.lh.1D.roi'];
end
sroi2v('brainmask+tlrc',roil,sid,'lh','sroi2v_lh_tlrc',suma,'Hua_SurfVol_ns_Alnd_Exp+tlrc');
for i=1:length(roir)
    roir{i} = [roir{i},'.rh.1D.roi'];
end
sroi2v('brainmask+tlrc',roir,sid,'rh','sroi2v_rh_tlrc',suma,'Hua_SurfVol_ns_Alnd_Exp+tlrc');

maskcalc('sroi2v_lh+orig[0]','step(a-2)',1,'vrois_lh');
maskcalc('sroi2v_rh+orig[0]','step(a-2)',1,'vrois_rh');

sroi2v('vrois_lh+orig','ecc_gt180.lh.1D.roi',sid,'lh','roi_eccgt180_lh',suma,'Hua_SurfVol_ns_Alnd_Exp+tlrc');
sroi2v('vrois_rh+orig','ecc_gt180.rh.1D.roi',sid,'rh','roi_eccgt180_rh',suma,'Hua_SurfVol_ns_Alnd_Exp+tlrc');