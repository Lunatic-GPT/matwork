function res=readInplaneRotAngle_uih(fname)

if isa(fname,'char')
    d=parseXML(fname);
else
    d=fname;
end
%%
root=get_child(d,'Root');
seq=get_name(root,'Seq');

gli=get_child(seq,'GLI');
sg=get_name(gli,'SliceGroup');
an=get_child(sg.Children(2),'InplaneRotAngleDisplay');

res=str2double(an(2).Children.Data);


function res=get_name(d,name)
res=[];
for i=1:length(d)

    if strcmp(d(i).Name,name)
        res=d(i);
        return;
    end

end

function res=get_child(d,name)
res=[];
for i=1:length(d.Children)

    if strcmp(d.Children(i).Name,name)
        res=d.Children(i).Children;
        return;
    end
end

function res=get_value(a)

if length(a)>1
    n=(length(a)-1)/2;
    res=zeros(1,n);
    for i=1:n
        res(i)=str2double(a(2*i).Children.Data);
    end
else
    n=(length(a.Children)-1)/2;
    res=zeros(1,n);
    for i=1:n
        res(i)=str2double(a.Children(i*2).Children.Data);
        
    end

    if n==1&&isnan(res)
      res=a.Children(2).Children.Data;
    end
end


