function myErrorBarCells(x,y,par,show_fit)
% y is a cell array which can be converted to mean and std

if ~exist('par','var')
    par=[];
end

if ~exist('show_fit','var')
    show_fit=false;
end


v=cell2meanarray(y);

ev=cell2stdarray(y);


h=myErrorBar(x,v,ev,par);


if show_fit
   
    
x2=repmat_match_cell2array(x,y);
v2=cell2array(y);

xx=[ones(length(x2),1),x2(:)];
b=xx\v2(:);

xr=min_max(x(:))'.*[0.5,2];
xx=[ones(length(xr),1),xr(:)];

yy=xx*b;
hold on;
plot(xx,yy,'r-','LineWidth',1);

xlim(min_max(x(:))'.*[0.9,1.1]);

 disp_corr(x2',v2);
 
 disp_corr(x2,v2,'Spearman');
end

set_plot(h,par);


