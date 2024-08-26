changebPar('reco','RECO_transposition',true,ones(1,24));
reco_rotate=readbPar('reco','RECO_rotate');
reco_rotate(2:end,1)=reco_rotate(1,1);
reco_rotate(2:end,2)=reco_rotate(1,2);
reco_rotate(2:end,3)=reco_rotate(1,3);
changebPar('reco','RECO_rotate',true,reco_rotate);


