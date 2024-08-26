function [a,ea]=read_results(fname,nfid)
%  [a,ea]=read_results(fname,nfid)
fid=fopen(fname,'r');

str=textscan(fid,'%s','Delimiter','\n');
a=[];
ea=[];
 for i=1:length(str{1})
   
   if strcmp(str{1}{i},'Amplitudes (-)')       
       for j=1:nfid
        a(j,:)=str2num(str{1}{i+j});
       end
   end
    if strcmp(str{1}{i},'Standard deviation of Amplitudes (-)')       
        for j=1:nfid
         ea(j,:)=str2num(str{1}{i+j});
        end
   end
   
 end

fclose(fid);
       
    