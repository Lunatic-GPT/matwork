function aa=TSE_angles_T2W_TYPE01_1000ms(necho,T2,T1,tau,aexc,evolPar)

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

tau2=evolPar(3)/necho;
TC1=evolPar(4);
TC2=evolPar(6);

Issa=zeros(1,necho);

for i=1:necho
    if i<=evolPar(3)/tau2
      tmp=sin(aexc)*exp(-tau2*i/TC1);
      if tmp>Iss
        Issa(i)=tmp;
      else
        Issa(i)=Iss;  
      end
    else        
      Issa(i)=Iss*exp(-tau*2*(i-necho_const)/TC2);
    end
end

aa=TSE_angles_Arb_S(Issa,T2,T1,tau,aexc);


