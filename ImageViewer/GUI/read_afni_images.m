function d=read_afni_images(fname)

    d=BrikLoad(fname);
    if size(d,4)>1
        
       answer=inputdlg('Select subrik','subbrik',1,{'0'});
       sb=str2num(answer{1});
    
    d=d(:,:,:,sb+1);
    end
    
    
   