function val = readPar(fid_dir,field)
%val = readPar(fid_dir or parfile,field)



    
 if exist([fid_dir,'.fid'],'dir') 
    fid_dir = [fid_dir,'.fid'];
 elseif exist([fid_dir,'.img'],'dir') 
    fid_dir = [fid_dir,'.img'];
 end

 if ~isdir(fid_dir) && exist(fullfile(pwd,fid_dir),'file')
    fname=fid_dir;
 else
    
 if exist(fullfile(fid_dir,'procpar.orig'),'file')
  fname=fullfile(fid_dir,'procpar.orig');
 else
   fname=fullfile(fid_dir,'procpar');
 end
end

%a=textscan(fid,'%s%s%s%s%s%s%s%s%s%s%s','delimiter',' ');
%a=textread(parfile,'%s%s%s%s%s%s%s%s%s%s%s',);
a=textread(fname,'%s');

ind = strmatch(field,a,'exact');
if length(ind) ~=1
    val=[];
    warning('Field not found');
    return;
end

type = a{ind+1};

ind = ind+11;
na = str2num(a{ind});
if isempty(na)
    error('Number of values error');
end

val = a{ind+1};

         if ~isempty(str2num(val))
            val = zeros(1,na);
              for i=1:na
               val(i) = str2num(a{ind+i});
              end
              
              if type == '6'
                  val = val/1000000;
              end
         elseif na==1
             val = a{ind+1};
         else
             val=a(ind+1:ind+na)';
         end
            

