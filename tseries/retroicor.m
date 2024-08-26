function retroicor(fname,resp,ecg,order)
 % retroicor(fname,resp,ecg,order)
 % the sampling rates are assumed to be 40 Hz for both resp and ecg.
    prefix = strtok(fname,'+');
    cmd = sprintf('3dretroicor -prefix pc%s ',prefix);
    if ~isempty(resp)
        cmd = sprintf('%s -resp %s -respphase respphase_%s.1D',cmd,resp,prefix);
    end
    
    if ~isempty(ecg)
        data = load(ecg);
        threshold = min(data)+(max(data)-min(data))*3/4;
        fprintf('ECG Threshold %4.3f\n',threshold);
        cmd = sprintf('%s -card %s -threshold %f -cardphase cardphase_%s.1D',cmd,ecg,threshold,prefix);
    end
    cmd = sprintf('%s %s',cmd,fname);
    unix(cmd);
    
    if ~isempty(ecg)
     phi_ecg = load(sprintf('cardphase_%s.1D',prefix));
    end
    
    if ~isempty(resp)
     phi_resp = load(sprintf('respphase_%s.1D',prefix));
    end
    
    [data,info] = BrikLoad(fname);
    nTR = size(data,4);
    TR = info.TAXIS_FLOATS(2);
    
    for k=1:size(data,3)
                t = TR*(0:nTR-1)+info.TAXIS_OFFSETS(k);
                if ~isempty(ecg)
                    p1 = spline(0.025*(0:length(phi_ecg)-1),phi_ecg,t);
                    p1 = p1(:);
                else
                    p1 = [];
                end
                if ~isempty(resp)
                    p2 = spline(0.025*(0:length(phi_resp)-1),phi_resp,t);
                    p2 = p2(:);
                else
                    p2 = [];
                end
                
                 nuas{k} = get_model(p1,p2,order);
    end
    
    for i=1:size(data,1)
        for j=1:size(data,2)
            for k=1:size(data,3)
                 data(i,j,k,:) = regress_out(data(i,j,k,:),nuas{k});
            end
        end
    end
    
    history = sprintf('retroicor(%s,%s,%s,%d)',fname,resp,ecg,order);
    WriteBrikEZ(data,info,history,['pc',prefix]);
    
         
    
    
function model = get_model(p1,p2,order)
        
model = [];
if ~isempty(p1)
    for i=1:order
      model = [model cos(i*p1) sin(i*p1)];
    end
end

if ~isempty(p2)
    
    for i=1:order
      model = [model cos(i*p2) sin(i*p2)];
    end
end


    
function ts = regress_out(ts,nuas)
        
        ts = ts(:);
 % model detrend       
        x = nuas - repmat(mean(nuas,1),[size(nuas,1),1]);
       
        b = lscov(x,ts);
        
        ts = ts - x*b;
        
%% this function is not used.  use the phase calculated from 3dretroicor        
function phi_ecg = get_phase_ecg(ecg)

    a = load(ecg);
    threshold = min(a)+(max(a)-min(a))*3/4;
    inda = find(a>threshold); %index for a
    
    inda_diff = diff(inda); %index difference of inda 

    ind_inda_diff = find(inda_diff>10);  % index for index difference.
    
    npk = length(ind_inda_diff)+1;
    
    pk_ind = zeros(1,npk);
    
    k = ind_inda_diff(1);
    [tmp,pk_ind(1)] = max(a(1:inda(k)));
    
    for j=1:npk
        if j==1
            k1 = 1;
        else
           k1 = ind_inda_diff(j-1);
        end
        if j==npk
            k2 =length(inda); 
        else
          k2 = ind_inda_diff(j);
        end
         [tmp,pk_ind(j)] = max(a(inda(k1):inda(k2)));
         pk_ind(j) = pk_ind(j)+inda(k1)-1;
    end
    
    phi_ecg = zeros(size(a));
    for j=1:npk
        if j==1
           t = (pk_ind(2)-pk_ind(1));
           phi_ecg(1:pk_ind(j)) = 2*pi*((1:pk_ind(1))-(pk_ind(1)-t))/t;
        else
            t = (pk_ind(j)-pk_ind(j-1));
            phi_ecg(pk_ind(j-1)+1:pk_ind(j)) = 2*pi*((pk_ind(j-1)+1:pk_ind(j))-pk_ind(j-1))/t; 
        end
    end
    
    if pk_ind(end)<length(a)  %  data after the last peak
            t = (pk_ind(end)-pk_ind(end-1));
            phi_ecg(pk_ind(end)+1:end) = 2*pi*((pk_ind(end)+1:length(a))-pk_ind(end))/t; 
    end
