function license_from_time(str,year,mon,day)

a(1:4)=digits_in_number(year,4);
a(5:6)=digits_in_number(mon,2);
a(7:8)=digits_in_number(day,2);
a=[a,a];

str=str+a;
fprintf('%c',str);
disp('');






function a=digits_in_number(n,nd)


for i=0:nd-1
a(i+1)=floor(n/10^(nd-i-1));
n=n-a(i+1)*10^(nd-i-1);
  
end
  