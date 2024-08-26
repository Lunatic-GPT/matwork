function scale = readCoilScalingFromProtocol(fname,nch)

 fid=fopen(fname,'r');
 scale=[];
 while ~feof(fid)
     d=fgetl(fid);
     
     if isempty(d)
         continue;
     end
     d=textscan(d,'%s');
     
     if nch==32
     coil_names = {'"A01"','"A02"','"A03"','"A04"','"A05"','"A06"','"A07"','"A08"',...
                   '"A09"','"A10"','"A11"','"A12"','"A13"','"A14"','"A15"','"A16"',...
                   '"A17"','"A18"','"A19"','"A20"','"A21"','"A22"','"A23"','"A24"',...
                   '"A25"','"A26"','"A27"','"A28"','"A29"','"A30"','"A31"','"A32"'};
     elseif nch==14
        coil_names = { '"H3P"', '"H4P"', '"H4S"','"H4T"','"H3S"','"H3T"','"NE2"',...
                       '"N2S"','"H1P"','"H2P"','"H2S"','"H2T"','"H1S"','"H1T"'};
     end
                   
     match=true;
     for i=1:nch
         
         if isempty(strmatch(coil_names{i},d{1}))
             match=false;
             break;
         end
         
     end
     
     if match
         
         scale=zeros(1,nch);
         for i=1:nch
             ind=strmatch(coil_names{i},d{1});
             scale(i)=str2double(d{1}{ind+6})+1i*str2double(d{1}{ind+9});
         end
         break;
         
     end
     
 end

 
 fclose(fid);