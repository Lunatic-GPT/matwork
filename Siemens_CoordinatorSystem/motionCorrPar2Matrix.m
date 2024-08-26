function motionCorrPar2Matrix
%%
% [cos(an),sin(an);-sin(an),cos(an)]
%% when transverse and rotation = 0 (readout: LR); 
% The comment par 1, 2, 3 = pFBData->dTransLog_0, pFBData->dTransLog_1,
% pFBData->dTransLog_2
% comment par 4, angle (in deg) for [1,2] as defined above (rotate around AP?)
% comment par 5, -angle for [0,2]
% comment par 6, angle for [0,1]
% Mat par 0, 1, 2 = A->P, L->R, I->S (trans)
%  


% when transverse and rotation = 90 (readout: AP)
% the comment par 1, 2, 3 = -pFBData->dTransLog_1, pFBData->dTransLog_0,
% pFBData->dTransLog_2 = R->L, A->P, I->S(need to check?)
% par 4, -angle for [1,2]
% par 5, -angle for [0,2]
% par 6, angle for [0,1]
% Mat par 0,1,2 = L->R, P->A, I-S


% dTransLog_0: phase
% dTransLog_1: readout
% dTransLog_2: par

% ImageComment (1): A->P (the amount of shift wrt the reference)
%              (2): L->R 
%              (3)  I->S


% issues:
% I-S orientation wrong
% L-R orientation correct
% A-P orientation also wrong

