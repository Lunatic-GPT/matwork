function co=var_conv(ts,n)

co= zeros(1,n);
for i=1:n
    m=length(ts);
     cov = corrcoef(ts(1:m-i+1),ts(i:m));
    % co(i)=mean(ts(1:m-i+1).*ts(i:m));
    co(i) = cov(1,2);
end
    