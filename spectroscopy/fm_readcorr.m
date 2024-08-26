function [m,cc]=fm_readcorr

if exist('/home/xiaopeng/labhome/vnmrsys/fastmap','dir')
    
file = '/home/xiaopeng/labhome/vnmrsys/fastmap/correlation';
else
    file='/data/users/xiaopeng/vnmrsys/fastmap/correlation';
end

data=textread(file,'%s');

ind=strmatch('reading',data,'exact');

if length(ind)~=8
    error('correlation file format error');
end

m=zeros(8,8);
cc=zeros(8,8);%component*shim coils
name={'x1','y1','z1','xz','yz','z2','xy','x2y2'};
i=1;
while i<length(data)
  if strcmp(data{i},'reading')
      i=i+6;
  end
  [coil,comp]=strtok(data{i},'-');
  ind = strmatch(coil,name,'exact');
  ind2 = strmatch(comp(3:end),name,'exact');
  if isempty(ind)
      error('no match');
  end
  
  m(ind2,ind)=str2double(data{i+3});
  cc(ind2,ind)=str2double(data{i+6});
  
  i=i+7;
      
end

