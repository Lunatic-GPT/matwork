function res=remove_mean_fullk(data,m)
% res=remove_mean_fullk(data,m)

t_dim=ndims(data);
      data=data.*m;
        tmp=sum(data,t_dim);
  
        ntmp=sum(m,t_dim);
        ntmp(ntmp==0)=1;
        mn=tmp./ntmp;
  
        sz=size(data);
        sz(1:end-1)=1;
        res=data-repmat(mn,sz);
        res=res.*m;
        