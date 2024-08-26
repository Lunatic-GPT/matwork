function dcm2dir_useHeader(root,plps)

% root: data will be saved in the root dir
% plps: phillips dicom, image type read from a private field (mag or
% phase).
warning off;
if ~exist('plps','var')
    plps=false;
end

if ~exist('root','var')
    root=pwd;
end


a=dir('*');

for i=1:length(a)
    
    
    if (a(i).isdir==1) && a(i).name(1)~='.'
        disp(a(i).name);
        cd(a(i).name);
        dcm2dir_useHeader(a(i).name,plps);
        cd ..;
        
    else
        
        try
            in=dicominfo(a(i).name);
            
            
            s=in.SeriesNumber;
            st=in.SeriesDescription;
            st=strrep(st,' ','_');
            %ind=find(st=='*');
            st(st=='*')=[];
            if isfield(in,'EchoNumber')
               echo=in.EchoNumber;
            else
                echo=1;
            end
            type='normalized';
            if plps
                type = in.Private_2005_140b;
            end
            
            if echo==1
                if strcmp(type,'normalized')
                    dname=[num2str(s),'_',st];
                else
                    dname=[num2str(s),'_',st,'_ph'];
                end
            else
                if strcmp(type,'normalized')
                    dname=[num2str(s),'_',st,'_echo',num2str(echo)];
                else
                    dname=[num2str(s),'_',st,'_echo',num2str(echo),'_ph'];
                end
            end
            
            if ~strcmp(filename(root),dname)  %not yet sorted             
                if ~exist(fullfile(root,dname),'dir')
                    mkdir(fullfile(root,dname));
                end
                movefile(a(i).name,fullfile(root,dname));
            end
            
        catch
            continue;
        end
    end
    
end

%cd(cur_dir);

