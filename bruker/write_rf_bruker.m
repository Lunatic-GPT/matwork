function write_rf_bruker(rfname,B1)
%[rf,tp]=read_rf_varian(rfname,B1max);



        a=fopen(rfname,'w');
 
      %  c=textscan(a,'%f%f','Delimiter',',','CommentStyle','##');
      for i=1:length(B1)
       fprintf(a,'%f,   %f\n',abs(B1(i)),mod(angle(B1(i))*180/pi,360));
      end
      fprintf(a,'##END=\n');
       fclose(a);
       
 
 