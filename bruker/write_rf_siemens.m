function write_rf_siemens(rf,rfname)
%[rf,tp]=read_rf_varian(rfname,B1max);



rf=rf/max(abs(rf));
        a=fopen(rfname,'w');
        
        for i=1:length(rf)
         fprintf(a,'%f  %f ; (%d)\n',abs(rf(i)),angle(rf(i)),i-1);
        end
        
        fclose(a);
        