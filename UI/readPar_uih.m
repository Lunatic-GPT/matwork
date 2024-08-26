function res=readPar_uih(fname,pname)

if isa(fname,'char')
d=parseXML(fname);
else
    d=fname;
end
%%
root=get_child(d,'Root');
seq=get_name(root,'Seq');
app=get_child(seq,'App');
a=get_name(app,pname);
n=(length(a.Children)-1)/2;
res=zeros(1,n);
for i=1:n
res(i)=str2double(a.Children(i*2).Children.Data);
end

function res=get_name(d,name)

for i=1:length(d)

    if strcmp(d(i).Name,name)
      res=d(i);
      return;
    end

end

function res=get_child(d,name)

for i=1:length(d.Children)

    if strcmp(d.Children(i).Name,name)
      res=d.Children(i).Children;
      return;
    end

end