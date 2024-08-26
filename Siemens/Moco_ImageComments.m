function res=Moco_ImageComments(dname,index_ref,index_plot)
% res=Moco_ImageComments(dname,index_ref,index_plot)
% res: 1 - 6: {'Roll (IS)','Pitch (RL)','yaw (AP)','IS','RL','AP'};
% same as dfile from dvolreg -dfile dfile
% index_ref: the reference time point (will be set as 0)
% index_plot: the time points to plot
a=readdPar(dname,'ImageComments',true);

mpar=zeros(length(a),6);
for i=2:length(a)
    
    b=textscan(a{i},'Motion: %f,%f,%f,%f,%f,%f');
    
    mpar(i,:)=cell2array(b);
    
end

mpar(:,[1,3,5])=-mpar(:,[1,3,5]);

[ori,rot]=get_orient_angle(dname);
%afni dfile: 1 - 6: {'Roll (IS)','Pitch (RL)','yaw (AP)','IS','RL','AP'};
% tested in MocoPar_vs_ImageComment

if ori=='T' && rot==0
% Tra; inplane 0 deg
% imgeComment: % 'PA', 'RL', 'SI', yaw(AP)', Pitch (LR), Roll (IS)'
% 
res=mpar(:,[6,5,4,3,2,1]);
elseif ori=='T' && int16(rot)==90

 % imgeComment: % 'PA', 'RL', 'SI', yaw(AP)', Pitch (LR), Roll (IS)'
  % res=mpar(:,[6,5,4,3,2,1]);
 %  res=mpar(:,[4,6,5,1,3,2]);
 %  res(:,2)=-res(:,2);
   res=mpar(:,[6,4,5,3,1,2]);
   res(:,[3,6])=-res(:,[3,6]);
elseif ori=='S' && int16(rot)==90
    res=mpar(:,[4,6,5,1,3,2]);
    res(:,[1,3,4,6])=-res(:,[1,3,4,6]);
elseif ori=='S' && int16(rot)==0
    res=mpar(:,[5,6,4,2,3,1]);
    res(:,[1,4])=-res(:,[1,4]);
    %res(:,[1,3,4,6])=-res(:,[1,3,4,6]);
elseif ori=='C' && int16(rot)==0
      res=mpar(:,[5,4,6,2,1,3]);  
elseif ori=='C' && int16(rot)==90
      res=mpar(:,[4,5,6,1,2,3]);  
      res(:,[2,5])=-res(:,[2,5]);
end


if nargout==0
    
    if exist('index_ref','var')
        res=res-res(index_ref,:);
    end
    
if ~exist('index_plot','var')
    index_plot=1:size(res,1);
end

plot_motion_afni(res(index_plot,:));

%saveas(gcf,sprintf('%s_moco',dname));
savetiffc(sprintf('%s_%dTo%d_moco',dname,index_plot(1),index_plot(end)));

end


function [ori,rot]=get_orient_angle(dcm)

%rot: deg
extp(dcm);
a(3)=readsPar([dcm,'.pro'],'sNavigatorArray.asElm[0].sCuboid.sNormal.dTra');
a(2)=readsPar([dcm,'.pro'],'sNavigatorArray.asElm[0].sCuboid.sNormal.dCor');
a(1)=readsPar([dcm,'.pro'],'sNavigatorArray.asElm[0].sCuboid.sNormal.dSag');
rot=readsPar([dcm,'.pro'],'sNavigatorArray.asElm[0].sCuboid.dInPlaneRot');

rot=180*rot/pi;
%copied from test_MOCO_Rotmat
%{ 
%% when rot =0; consistent with sROT_MATRIX(orient,0); %2/25/2021
if a(3)==1 %transverse
    Mf0=[0,-1,0;...
        -1,0,0;...
        0,0,-1];  %PE, RO, PAR
elseif a(2)==1 %coronal
    Mf0=[1,0,0;...
        0,0,-1;...
        0,1,0];
else % sagittal
    Mf0=[0,0,1;...
        -1,0,0;...
        0,-1,0];
end
%}

[~,ind]=max(abs(a));

if ind==3 %transverse
    ori='T';  
elseif ind==2 %coronal
    ori='C';
else % sagittal
   ori='S';
end