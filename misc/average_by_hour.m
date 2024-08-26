function [t,out,eout]=average_by_hour(time,step,data)
%[t,out,eout]=average_by_hour(time,step,data)
% time and data are cell arrays
% each data element should have the same number of rows as the time.
% step is the time step for combining the data within each step.

time_max=0;
nd=0;
for i=1:length(time)
    if time_max<max(time{i})
        time_max=max(time{i});
    end
    
    if size(data{i},2)>0
        nd=size(data{i},2);
    end
end

nh=ceil(time_max/step);



out=zeros(nh,nd);
eout=zeros(nh,nd);
t=((1:nh)-0.5)*step;

sdata=cell(1,nh);
for i=1:nh
    tm=[];
  %  etm=[];
    for j=1:length(time)
        
        ind=find(time{j}<=i*step&time{j}>(i-1)*step);
        if ~isempty(ind)
          tm=[tm;data{j}(ind,:)];
   %       etm=[tm,data{j}(ind)];
        end
        
        
    end
    
    if length(tm)>0
      out(i,:)=mean(tm,1);
      eout(i,:)=std(tm,0,1)/sqrt(size(tm,1));
    end
    
    sdata{i}=tm;
        
end

out(out==0)=NaN;
eout(isnan(out))=NaN;

for i_start=1:nh
 if ~isempty(sdata{i_start})
     break
 end
end


for i=i_start+1:nh
    x1=sdata{i_start};
    x2=sdata{i};
    x1(x1==0)=NaN;
    x2(x2==0)=NaN;
    if ~any(~isnan(x1)) || ~any(~isnan(x2)) 
        continue;
    end
      [h,p]=ttest2(x1,x2);
      fprintf('%d - %d: p=%4.3f\n',i_start,i,p);
    
    
end

    
    
    
    
    

