function pos=read_swc(fname,dm)
% read trace files from simple neurite tracer into an list, one element for
% each trace
% dm: dimension in mm

if exist(fname,'dir')
    dir_str=dir(fullfile(fname,'*.swc'));
    pos=cell(1,length(dir_str));
    for j=1:length(dir_str)
        fname=fullfile(dir_str(j).folder,dir_str(j).name);
        fid=fopen(fname,'r');
              
        a=textscan(fid,'%f %f %f %f %f %f %f','CommentStyle','#');
        fclose(fid);
        pos{j}=zeros(length(a{1}),3);
        for i=1:length(a{1})
            pos{j}(i,:)=[a{3}(i),a{4}(i),a{5}(i)]./dm+1;
        end
    end
    
else
     fid=fopen(fname,'r');        
     a=textscan(fid,'%f %f %f %f %f %f %f','CommentStyle','#');
     fclose(fid);
     j=0;
     ind=0;
     for i=1:length(a{1})
     
        if a{7}(i)<0
            j=j+1;
            ind=1;
        end
        pos{j}(ind,:)=[a{3}(i),a{4}(i),a{5}(i)]./dm+1;    
    
        ind=ind+1;
     end
end
