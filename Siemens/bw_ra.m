function [bw,refa]=bw_ra(dname)
extp(dname);
bw=readdPar(dname,'PatientWeight');
refa=readsPar([dname,'.pro'],'flReferenceAmplitude');
%fprintf('Weight = %3.0f; ReferenceAmplitude =  %3.0f\n',bw,refa);

function extp(fname)
%% extract protocol

if exist(fname,'dir')
dir_str=dir(fname);

i=0;
while 1
    if exist(fullfile(dir_str(end-i).folder,dir_str(end-i).name),'dir')
        
        i=i+1;
        
    else
        fid=fopen(fullfile(dir_str(end-i).folder,dir_str(end-i).name),'rb');
        break;
    end
end
prefix=strtok(fname,'.');
else
fid=fopen(fname,'rb');

prefix=strtok(fname,'.');

end
%for i=1:length(s)
    
 %   if strmatch(s{i},'ASCCONV')
 fid2=fopen([prefix,'.pro'],'w');
 start=0;
while ~feof(fid)
    s=fgetl(fid);
    if ~isempty(strfind(s,'ASCCONV'))
        start=start+1;
        
        if start==1
            fprintf(fid2,'### ASCCONV BEGIN ###\n');
            continue;
        else
            fprintf(fid2,'### ASCCONV END ###\n');
            break;
        end
        
    end
    
    if start==1
        
        fprintf(fid2,'%s\n',s);
        
    end
end

fclose(fid);
fclose(fid2);


        

function res3=readdPar(dname,par,allfile)
%res=readdPar(dname,par,allfile)

d=dir2(dname);

if ~exist('allfile','var')  || ~allfile
    
    for i=1:length(d)
        if isdir(fullfile(d(i).folder,d(i).name))
            continue;
        end
        in=dicominfo(fullfile(d(i).folder,d(i).name));
        
        break;
    end
    if exist('par','var')
        res3=getfield(in,par);
    else
        
        disp(in);
    end

else
    
    res3={};
    for i=1:length(d)
     %   disp(i);
        in=dicominfo(fullfile(d(i).folder,d(i).name));
        try
            res3{i}=getfield(in,par);
        catch
            res3{i}='';
        end
        
    end
    
    
    if isa(res3{1},'char')
       return; 
    end
    try
     res3=cell2mat(res3);
    catch
        
    end
    
end

function res3=readsPar(fname,par)
%res=readbPar(fname,par,isnum)

fid=fopen(fname,'r');
res3=[];
while ~feof(fid)
    
  b=fgetl(fid);
 

  ind=strfind(b,par);
  if ~isempty(ind) 
    res=b(ind+length(par):end);
    res=strrm(res,' ');
    res=strrm(res,'=');

  else
      continue;
  end


res2=str2double(res);

if ~isnan(res2)    
      res3(end+1)=res2;
else
    res3{end+1}=res;
end
end
if isempty(res3) && nargout==0
    disp([par, ' not found']);
end

fclose(fid);

