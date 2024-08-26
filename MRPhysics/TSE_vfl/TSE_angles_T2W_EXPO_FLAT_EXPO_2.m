function [aa_final,Issa_final]=TSE_angles_T2W_EXPO_FLAT_EXPO_2(necho,T2,T1,esp,aexc,evolPar,faMax)

%
% adEvolutionParameters[sEvolution.ulName][6]: total time for the second
% exp
% adEvolutionParameters[sEvolution.ulName][5]: time const for the second exp
% adEvolutionParameters[sEvolution.ulName][4]: total time for the first exp
% adEvolutionParameters[sEvolution.ulName][3]: time const for the first exp

%{ 0.5, 100.0, 1000.0, 174.7030, 119.4969, 289.8299, 465.4088, 0.0, 0.0, 0.0 }

%
%fa=TSE_angles_exp_flat_exp(Iss,necho,T2,T1,tau,aexc,40*tau,40*tau,necho_const);

%evolPar=[0.5, 100.0, 1000.0, 174.7030, 119.4969, 289.8299, 465.4088];

aexc_rad=aexc*pi/180;
TC1=evolPar(6);
TC2=evolPar(6);

Issa=zeros(1,necho);
n1=evolPar(5);
n2 = floor(necho-evolPar(7)/esp);
b=evolPar(8); %exponent
%s0=sin(aexc_rad);
%% need to check the sequence diagram to be sure 4/8/2020

for s0=0.01:0.01:1
    for i=1:necho
        if i<=n1
            Issa(i)=s0*exp(-(esp*i)^b/TC1);
        elseif i<=n2
            Issa(i)=Issa(n1);
            
        elseif i>n2
            Issa(i)=Issa(n2)*exp(-esp*(i-n2)/TC2);
        end
    end
    
    aa=TSE_angles_Arb_S(Issa,T2,T1,esp/2,aexc);
    
    if length(aa)==length(Issa) && max(aa)<faMax
        aa_final=aa;
        Issa_final=Issa;
    else
        break;
    end
    
end
%figure(3);hold on;plot(aa,'r');

%disp('');
