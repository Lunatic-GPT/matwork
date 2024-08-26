function rvtc(fname,resp)
% rvtc(fname,resp)
    d = load(resp);
    [maxtab,mintab] = peakdet(d,300);
    
    fprintf('Resperatory rate %3.2f/min\n',size(maxtab,1)*60/(length(d)*0.025));
    rvt = [];
    t_rvt = [];
    
    if mintab(1,1) < maxtab(1,1)
     imaxtab = 1;
    else
     imaxtab = 2;
    end
    
    for i=1:size(mintab,1) - 1
    
      if mintab(i+1,1)<maxtab(imaxtab,1)
          error('two consecutive valleys without a peak');
      end
          
      rvt(end+1) = (maxtab(imaxtab,2)-mintab(i,2))/(mintab(i+1,1)-mintab(i,1));
      imaxtab = imaxtab + 1;
      t_rvt(end+1) = (mintab(i+1,1)+mintab(i,1))/2*0.025-0.025;
    end
    
     
    [data,info] = BrikLoad(fname);
    nTR = size(data,4);
    TR = info.TAXIS_FLOATS(2);
    
    model = cell(size(data,3),51);
    for k=1:size(data,3)
        for i = 1:51
         t = TR*(0:nTR-1)+info.TAXIS_OFFSETS(k)+i-11;
         model{k,i} = interp1(t_rvt,rvt,t,'linear');   
        end
    end
            
    delay = zeros(size(data,1),size(data,2),size(data,3));
    
    for i=1:size(data,1)
        for j=1:size(data,2)
            for k=1:size(data,3)
               [data(i,j,k,:),delay(i,j,k)] = regress_out(data(i,j,k,:),model(k,:)); 
            end
        end
    end
    
    prefix = strtok(fname,'+');
    history = sprintf('rvtc(%s,%s)',fname,resp);
    WriteBrikEZ(data,info,history,['rvtc_',prefix]);
    WriteBrikEZ(delay,info,history,['rvtcDelay_',prefix]);
    
    
function [ts,delay] = regress_out(ts,nuas)
        
        ts = ts(:);
 % model detrend
        var_fit = zeros(1,length(nuas));
        
        for i=1:length(nuas)
          ind = ~isnan(nuas{i});
          ts_tmp = ts(ind);
          x_tmp = nuas{i}(ind);
          x_tmp = x_tmp - mean(x_tmp);
          
          b = lscov(x_tmp(:),ts_tmp(:));
          var_fit(i) = sum((x_tmp*b).^2);
          
        end
        
        [tmp,ind_max] = max(var_fit);
        ind = ~isnan(nuas{ind_max});
          ts_tmp = ts(ind);
          x_tmp = nuas{ind_max}(ind);
          x_tmp = x_tmp - mean(x_tmp);
          x_tmp = x_tmp(:);
          b = lscov(x_tmp,ts_tmp(:));
        ts(ind) = ts(ind) - x_tmp*b;
        delay = ind_max-11;
        
        
        
        
        
    