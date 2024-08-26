function m=read_mat_roi(prefix,ns)

       
       for i=1:ns
           
           
       fname=sprintf('%sX%d.mat',prefix,i);
       
        if ~exist(fname,'file') 
            if exist('m','var')
                m(:,:,i)=0;
            end
           continue;
        else
           tmp=load(fname);    
           m(:,:,i)=tmp.roi;
        end
       
       end
       
       
       
       
