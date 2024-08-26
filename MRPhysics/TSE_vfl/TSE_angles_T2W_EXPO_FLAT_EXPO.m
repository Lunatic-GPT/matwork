function aa=TSE_angles_T2W_EXPO_FLAT_EXPO(necho,T2,T1,esp,aexc,evolPar)

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

TC1=evolPar(4);
TC2=evolPar(6);

Issa=zeros(1,necho);
n1=floor(evolPar(5)/esp);
n2 = floor(necho-evolPar(7)/esp);
% need to check the sequence diagram
for i=1:necho
    if i<=n1
      Issa(i)=sin(aexc)*exp(-esp*i/TC1);
    elseif i<=n2
      Issa(i)=Issa(n1);
        
    elseif i>n2 
      Issa(i)=Issa(n2)*exp(-esp*(i-n2)/TC2);
    end
end

aa=TSE_angles_Arb_S(Issa,T2,T1,esp/2,aexc);


