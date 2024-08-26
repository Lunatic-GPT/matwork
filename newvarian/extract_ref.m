function ref=extract_ref(fname,i,prefix)

% extract_ref(fname,i,prefix)
% i is a 2*3 matrix specifying the coordinates of the upper left and lower
% right corners, 0 based;
if isa(fname,'char')
    
  data=rdSdt(fname);
else
    data=fname;
end

    ind1 = i(1,1):i(2,1);
    ind2 = i(1,2):i(2,2);
    ind3 = i(1,3):i(2,3);
    tmp = data(ind1+1,ind2+1,ind3+1,:);
    tmp = mean(tmp,1);
    tmp = mean(tmp,2);
    tmp = mean(tmp,3);
    
    ref=squeeze(tmp);
    figure;plot(ref);
    
    save(prefix,'ref','-ascii');
    