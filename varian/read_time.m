function t = read_time(fid_prefix,start_end)
%t = read_time(fid_prefix,start_end)
% t is a 2*3 vector. containing the start and end times in [hr,min,sec].
% start_end: 1 to get only the start time, 2 to get only the end time.
% 3 to get the average. 0 to get both (default)

if ~exist('start_end','var')
    start_end=0;
end

fid = fopen([fid_prefix,'.fid/log']);
str=textscan(fid,'%s%s%s%s%s%s%s%s');

t=zeros(2,3);

    
 t(1,1)=str2double(str{4}{1}(1:2));
 t(1,2)=str2double(str{4}{1}(4:5));
 t(1,3)=str2double(str{4}{1}(7:8));
 

 t(2,1)=str2double(str{4}{end}(1:2));
 t(2,2)=str2double(str{4}{end}(4:5));
 t(2,3)=str2double(str{4}{end}(7:8));
 

fclose(fid);

if start_end==1
    t=t(1,:);
elseif start_end==2
    t=t(2,:);
elseif start_end==3
    t=add_time(t(1,:),(t(2,:)-t(1,:))/2);
end
