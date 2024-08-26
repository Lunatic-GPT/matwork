function [total,an,shift]=afni_motionPar_newBase(dfile,ref)
% ref: the reference time point; 1 based
% results are in degrees
% newBase from afni_motionPar_newBase

m_all=afni_motionPar2Mat(dfile);
nt=size(m_all,1);

m_all=reshape(m_all,[nt,4,3]);
m_all=permute(m_all,[3,2,1]);
m_all(end+1,:,:)=repmat([0,0,0,1],[1,1,nt]);

m_ref=m_all(:,:,ref);
for i=1:nt

    m_all(:,:,i)=m_ref\m_all(:,:,i);
    total(:,i)=get_total(m_all(1:3,1:3,i));    
    an(:,i)=get_rot_axis(m_all(1:3,1:3,i));
    
end

shift=m_all(1:3,4,:);

function res=get_total(m)

fitfunc=@(x) vec(transform_matrix_rotation_arb_axis(x(1),x(2),x(3))-m);
a=optimoptions('lsqnonlin','Display','none','OptimalityTolerance',0,'FunctionTolerance',0,'StepTolerance',1e-8);
res=lsqnonlin(fitfunc,[0,0,0],[],[],a);



function res=get_rot_axis(m)

fitfunc=@(x) vec(mat_rot_axis(x)-m);
a=optimoptions('lsqnonlin','Display','none');
res=lsqnonlin(fitfunc,[0,0,0],[],[],a);



function m=mat_rot_axis(an)
    m_IS = transform_matrix_rotation_arb_axis(0,0,an(1));
    m_RL = transform_matrix_rotation_arb_axis(90,0,an(2));
    m_AP = transform_matrix_rotation_arb_axis(90,90,an(3));
    
    %m=m_IS*m_RL*m_AP;
    %m=m_IS*m_AP*m_RL;
    %m=m_AP*m_RL*m_IS;
    %m=m_AP*m_IS*m_RL;
    %m=m_RL*m_IS*m_AP;
    m=m_AP*m_RL*m_IS;  % not sure which one is correct
    
    
    

