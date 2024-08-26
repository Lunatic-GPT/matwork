function ts_detrend_cstm(dname,mask,prefix)
% ts_detrend(dname,mask,prefix)
% ts_detrend: remove the constant and linear terms in the single trial
% average time
% tps: the time points used to estimate the baseline. one based
% order: the polynomial order to fit the base line.
% prefix: default [prefix of data, '_dtr'];

tic;

fprintf('detrending %s ...\n',dname);
        [data,info] = BrikLoad(dname);
        sz = size(data);
        [nv,bl] = nv_roi(dname,mask);
        d2 = zeros(size(data));
        for i=1:size(data,1)
            for j=1:size(data,2)
                for k=1:size(data,3)
                    ts = squeeze(data(i,j,k,:));
                    
                    d2(i,j,k,:) = ts_detrend_cstm1D(ts(:),bl(:),-1);     
                end
            end
        end
     
        
        history = sprintf('\\nts_detrend_cstm(%s,mask = %s,prefix=%s)',dname,mask,prefix);
        info.HISTORY_NOTE = [info.HISTORY_NOTE,history];
        info.BRICK_TYPES = ones(1,sz(4))*3; %float type
        info.BRICK_FLOAT_FACS = zeros(1,sz(4))*3;
        opt.OverWrite = 'y';
        opt.Prefix = prefix;
        
        
        WriteBrik(d2,info,opt);
        
        disp([mfilename ' finish in ', num2str(toc), ' s']);
        