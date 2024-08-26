function [rf,tp]=read_rf_varian(rfname,pw,B1max)
%[rf,tp]=read_rf_varian(rfname,pw,B1max);

root=pwd;
    if isempty(rfname)
        a_tmp=repmat([0,0,1],[100,1]);     
    elseif strcmp(rfname,'x')
        a_tmp = [-pi/2,B1max,1];
    elseif strcmp(rfname,'y')
        a_tmp = [0,B1max,1];
    elseif strcmp(rfname,'-x')
        a_tmp = [pi/2,B1max,1];
    elseif strcmp(rfname,'-y')
        a_tmp = [pi,B1max,1];
    else
        a_tmp=textread(fullfile(root,rfname),'','commentstyle','shell');
        a_tmp(:,1)=a_tmp(:,1)*pi/180;
        a_tmp(:,2) = a_tmp(:,2)/max(a_tmp(:,2))*B1max;
    end
    
    
     t = a_tmp(:,3);
     tp = t/sum(t)*pw;
 rf = a_tmp(:,2).*(cos(a_tmp(:,1))+1i*sin(a_tmp(:,1)));