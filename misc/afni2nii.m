function afni2nii(dname)

[brikdata,info]=BrikLoad(dname);

brikdata=reshape(brikdata,[sz(1:2),prod(sz(3:4)),sz(5)]);
%save(dname,'brikdata');
     
%% save the average time series for each trial type.    
vsize=abs(info.DELTA);

writeanalyze(brikdata,dname,vsize,'int16');
    disp([mfilename ' finish in ', num2str(toc), ' s']);
