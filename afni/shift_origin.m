function shift_origin(fname,shift)
% shift_origin(fname,shift)
% shift should be a 1*3 vector along R(-)-L(+), A(-)-P(+); and I(-)-S(+).


[err,b]=BrikInfo(fname);

orient=b.ORIENT_SPECIFIC;
origin=b.ORIGIN;
ind=ceil((orient+1)/2);

disp(origin+shift(ind));



