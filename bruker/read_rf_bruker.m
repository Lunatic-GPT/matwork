function rf=read_rf_bruker(rfname,B1max)
%[rf,tp]=read_rf_varian(rfname,B1max);



        a=fopen(rfname);
        c=textscan(a,'%f%f','Delimiter',',','CommentStyle','##');
          fclose(a);
       
        a_tmp(:,2)=c{2}*pi/180;
        a_tmp(:,1) = c{1}/max(c{1})*B1max;
    
 rf = a_tmp(:,1).*(cos(a_tmp(:,2))+1i*sin(a_tmp(:,2)));
 
