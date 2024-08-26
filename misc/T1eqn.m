function Mfinal=T1eqn(Mint,Mss,T1,t)
%Mfinal=T1eqn(Mint,Mss,T1,t)
Mfinal=Mss+(Mint-Mss)*exp(-t/T1);