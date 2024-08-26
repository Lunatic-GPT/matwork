function res=plot_range(d,margin)

ymin=min(d);
ymax=max(d);


    dmin=-margin*(ymax-ymin);   
    dmax=margin*(ymax-ymin);


res=[ymin+dmin,ymax+dmax];


