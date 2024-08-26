function [factor]=get_corr(fid_water,fid_met,lsfid)
% factor=get_corr(fid_water,fid_met,lsfid)
% correction factor used for normalizing water-suppressed spectra by water signal. 
% s_new=s_old*factor, s_new is the ratio of metabolite to water
% signals under same gain and nt.

    fw=read_fid(fid_water);
    nt=readPar(fid_met,'nt');
    ntw=readPar(fid_water,'nt');
    g=readPar(fid_met,'gain');
    gw=readPar(fid_water,'gain');

    fw=fw-repmat(mean(fw(end-100:end,:),1),[size(fw,1),1]);


    factor=1/fw(lsfid+1)/nt*ntw/db2factor(g-gw);
    factor=abs(factor);

