function data=scale2int16(data)

mx=max(abs(data(:)));

data=int16(data/mx*(2^15-1));



