a=load('DspData_RFD.dsv');

cval=0;
i=1;
d=[];
while 1
    
    if a(i)==0 && i>1&&a(i-1)==0
        
        d(end+1:end+a(i+1)+1)=cval;
        i=i+2;
        
    else
        
        cval=cval+a(i);
        d(end+1)=cval;
        i=i+1;
    end
   
    
    if i>length(a)
        break;
    end
end

        