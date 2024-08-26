function version=rawDataVersion(filename)




fid = fopen(filename,'r','l','US-ASCII');

% start of actual measurement data (sans header)
fseek(fid,0,'bof');

firstInt  = fread(fid,1,'uint32');
secondInt = fread(fid,1,'uint32');

% lazy software version check (VB or VD?)
if and(firstInt < 10000, secondInt <= 64)
    version = 'vdve';
  
else
    % in VB versions, the first 4 bytes indicate the beginning of the
    % raw data part of the file
    version  = 'vb';
   
end

fclose(fid);


