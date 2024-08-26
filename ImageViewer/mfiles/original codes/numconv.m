function filename = numconv(fnum)
% convert number to four digits
% takes a number and returns a string
fnum_str = num2str(fnum);
if fnum < 10
    fnum_str = ['000' fnum_str];
elseif (9 < fnum) && (fnum < 100)
    fnum_str = ['00' fnum_str];
elseif (99 < fnum) && (fnum < 1000)
    fnum_str = ['0' fnum_str];
end
filename = fnum_str;