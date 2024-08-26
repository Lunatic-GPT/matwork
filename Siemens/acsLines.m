function a=acsLines(sk)


a=[];
for i=1:length(sk)-1
    if abs(sk(i)-sk(i+1))==1
        a(end+1)=sk(i);
    elseif i>1&&abs(sk(i)-sk(i-1))==1
        a(end+1)=sk(i);
    end
end
