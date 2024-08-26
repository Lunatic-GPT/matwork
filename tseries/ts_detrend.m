function ts = ts_detrend(dname,tps,order,prefix)
% ts_detrend(dname,tps,order,prefix)
% ts_detrend: remove the linear terms in the single trial
% average time
% dname: afni file name or a 2d/4d matrix. For a 2d matrix, time is the
% second dimension (row) and detrending is performed for each row.
% tps: the time points used to estimate the baseline. one based
% order: the polynomial order to fit the base line. 0,1,2...
% prefix: default [prefix of data, '_dtr']; not used if dname is a 2d matrix.

tic;


if ~exist('prefix','var')  && isa(dname,'char')
   prefix = strtok(dname,'+');
   prefix = [prefix,'_dtr'];
end

         if isa(dname,'char')
             fprintf('detrending %s ...\n',dname); 
             [data,info] = BrikLoad(dname);
         else
             data = dname;
         end

         nd2=false;
     if ndims(data)==2
         nd2=true;
         
         if size(data,2)==1
             
           data=reshape(data,[1,1,1,size(data,1)]);
         else
           data=reshape(data,[size(data,1),1,1,size(data,2)]);
         end
         
     end
     

       sz = size(data);
       if order >0  
        for i=1:size(data,1)
            disp(i);
            for j=1:size(data,2)
                for k=1:size(data,3)
                    
                    ts = squeeze(data(i,j,k,tps));
                    p = polyfit(tps,ts',order);
                    t = 1:sz(4);
                    
                    bl = zeros(1,sz(4));
                    for ip=0:order
                       bl = bl+t.^(order-ip)*p(ip+1);
                    end
                    bl = shiftdim(bl,-2);                       
                    data(i,j,k,t) = data(i,j,k,t) - bl+mean(data(i,j,k,t),4); 
                end
            end
         %   disp(i);
        end
       else 
          data = data - repmat(mean(data(:,:,:,tps),4),[1,1,1,sz(4)]);
       end
        ts = data;
        if exist('prefix','var')
            
          save(prefix,'ts');
          
        elseif nd2
            ts=squeeze(ts);
            
        end 
        
        
        
        