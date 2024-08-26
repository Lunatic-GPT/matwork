function res=xformMat2totalAngle(m)

fitfunc=@(x) vec(transform_matrix_rotation_arb_axis(x(1),x(2),x(3))-m);
a=optimoptions('lsqnonlin','Display','none','OptimalityTolerance',0,'FunctionTolerance',0,'StepTolerance',1e-8);
res=lsqnonlin(fitfunc,[0,0,0],[],[],a);


    
    
    

