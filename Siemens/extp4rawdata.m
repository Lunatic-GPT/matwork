function extp4rawdata(filename)
%[actualData,ushLine,ushPartition,ushSlice] = readMeasDat(filename,max_nlines,offset[,remove_os])
% filename: data file name
% max_nlines; the max number of lines to read
% offset: the starting offset; 0 start from the first line. 
 % remove_os: If true, remove the 2* oversampling along readout.  If remove_os is true, then actualData is already Fourier
 % transformed along readout. default: false.
 % ushLine and ushPartition are 0 based;
 

fid = fopen(filename,'r','l','US-ASCII');

% start of actual measurement data (sans header)
fseek(fid,0,'bof');

firstInt  = fread(fid,1,'uint32');
secondInt = fread(fid,1,'uint32');

% lazy software version check (VB or VD?)
if and(firstInt < 10000, secondInt <= 64)
    
    disp('Software version: VD/VE (!?)');

    % number of different scans in file stored in 2nd in
    % measOffset: points to beginning of header, usually at 10240 bytes
    measOffset = fread(fid,1,'uint64');
    fseek(fid,measOffset,'bof');
    hdrLength  = fread(fid,1,'uint32');

else
    % in VB versions, the first 4 bytes indicate the beginning of the
    % raw data part of the file
    disp('Software version: VB (!?)');
    hdrLength  = firstInt;
end


 xprotocol=[strtok(filename,'.'),'.xprotocol'];
  
 if exist(xprotocol,'file')
     fseek(fid, hdrLength, 0);
 else
     header=fread(fid,hdrLength,'*char');
     fid2=fopen(xprotocol,'w');
     
     for i=5:length(header)-12
         fprintf(fid2,'%c',header(i));
     end
     fclose(fid2);
     
     if ~exist([strtok(filename,'.'),'.pro'],'file')
         extp(xprotocol);
     end
 end



 fclose(fid);


