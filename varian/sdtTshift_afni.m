function sdtTshift_afni(fname,tpattern)

sdt2afni(fname);

tr = readPar(fname,'tr');

cmd= sprintf('3dTshift -tpattern %s -TR %f -prefix %s_tShft %s+orig',tpattern,tr,fname,fname);

unix(cmd);
afni2sdt([fname,'_tShft'],false);

