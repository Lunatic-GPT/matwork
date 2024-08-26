dir_str=dir('*.dat');

for i=1:length(dir_str)

    prefix=strtok(dir_str(i).name,'.');
    
    
    nois=readMeasDat([prefix,'.dat'],32,0,true);  

    save(prefix,'nois');
    
end