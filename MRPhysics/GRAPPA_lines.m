function res=GRAPPA_lines(acs,npe,pat)

res=1:pat:npe;
res=[res,npe/2-acs/2+1:npe/2+acs/2];
res=unique(res);




