function rf=read_rf_siemens(rfname)
%[rf,tp]=read_rf_varian(rfname,B1max);



        a=fopen(rfname);
        d=[];
        i=1;
        while 1
            c=fgetl(a);
            if c==-1
                break;
            end
         
            if isempty(c)
                continue;
            end
            
            ind=findstr(c,':');
            if isempty(ind) && length(c)>1
                c=strtok(c,';');
                
                d(i,:)=str2num(c);
                i=i+1;
            end
            
        end
        fclose(a);
        
 rf = d(:,1).*(cos(d(:,2))+1i*sin(d(:,2)));