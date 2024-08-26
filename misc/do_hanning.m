function res=do_hanning(data,dim)

n=size(data,dim);

flt=shiftdim(hanning(n),dim-1);

res=data.*flt;
