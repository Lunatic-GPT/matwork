 a=read_rf_bruker('c:/Users/xiaopeng/Dropbox/ForBruker/wave/sinc.exc',1);
 
fprintf('{');
for i=1:1024
    
    if mod(i,20)~=0
    fprintf('%4.3f,',abs(a(i)));
    else
        fprintf('%4.3f,\n',abs(a(i)));
    end
    
    
end


fprintf('}');

fprintf('{');
for i=1:1024
    
    
        if a(i)>0
           fprintf('0.0,');
        else
           fprintf('3.14159,');
        end
      if mod(i,20)==0
          fprintf('\n');
      end
        
    
end


fprintf('}');
