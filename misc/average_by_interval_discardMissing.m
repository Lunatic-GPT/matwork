function [t,out,eout,n,p,exclude,d,d2]=average_by_interval_discardMissing(time,intervals,data)
%[t,out,eout,n,p,exclude,d,d2]=average_by_interval_discardMissing(time,intervals,data)
% time and data are cell arrays
% each data element should have the same number of rows as the time.
% all data element should have the same number of columns.
% step is the time step for combining the data within each step.
% d2: subjects with missing data are not excluded.

for i=1:length(data)
    if ~isempty(data{i}) && size(data{i},1)==1
       data{i}=data{i}';  %convert to column vector
    end
end


for i=1:length(data)
    if ~isempty(data{i})
       nd=size(data{i},2);
       break;
    end
end

nh=length(intervals)-1;

ns=length(time);
d=nan(nh,nd,ns);
t=0.5*(intervals(1:end-1)+intervals(2:end));

for i=1:ns
    for j=1:nh
        
        ind=find(time{i}<=intervals(j+1)&time{i}>intervals(j));
        if ~isempty(ind)
          d(j,:,i)=mean(data{i}(ind,:),1);     
     
        end
        
        
    end
    
end

exclude=false(1,ns);
for i=1:ns
    tmp=d(:,:,i);
    if any(isnan(tmp(:)))
        exclude(i)=true;
    end
end

d2=d;    
d(:,:,exclude)=[];
out = mean(d,3);
eout=std(d,0,3)/sqrt(size(d,3));
n=size(d,3);
p=ones(1,nh+1);
if nd==1
    
    if size(d,3)>2
      for i=2:nh
        p(i)=ttest_1vec(d(1,1,:)-d(i,1,:));
   %    fprintf('1 - %d: p=%4.3f\n',i,p(i));  
    
      end
      
      
      p(nh+1)=ftest_matrx(squeeze(d(:,1,:))');
    end
end

    
    
    
    

