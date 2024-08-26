function y=do_rplc(y,y_rplc)


sel=~isnan(y_rplc);
y(sel)=y_rplc(sel);
