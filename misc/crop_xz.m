function res = crop_xz(x,s)


    m = size(x);
    if length(s) < length(m)
	    s = [s, m(length(s)+1:end)];
    end
	
    if sum(m==s)==length(m)
	res = x;
	return;
    end

    
    for n=1:length(s)
	    idx{n} = floor(m(n)/2)+1+ceil(-s(n)/2) : floor(m(n)/2)+ceil(s(n)/2);
    end

    % this is a dirty ugly trick
    cmd = 'res = x(idx{1}';
    for n=2:length(s)
    	cmd = sprintf('%s,idx{%d}',cmd,n);
    end
    cmd = sprintf('%s);',cmd);
    eval(cmd);





