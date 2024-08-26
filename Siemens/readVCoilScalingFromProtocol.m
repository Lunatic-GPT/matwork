function scale = readVCoilScalingFromProtocol(fname)

 fid=fopen(fname,'r');
 scale=[];
 while ~feof(fid)
 d=fgetl(fid);
 
 if isempty(d)
     continue;
 end
 d=textscan(d,'%s');
 
 if ~isempty( strmatch(sprintf('"VC"',i),d{1}))
     
        ind=strmatch('"VC"',d{1});
        scale=str2double(d{1}{ind+6})+1i*str2double(d{1}{ind+9});
         
    break;
 end
 
 
 end

 
 fclose(fid);