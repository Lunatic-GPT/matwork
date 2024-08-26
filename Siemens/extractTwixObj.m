function res=extractTwixObj(obj)
sz=obj.dataSize;
n=length(obj.memPos);
sz=[sz(1:2),n];

Data=obj.readData(obj.memPos,1:n,1:n,{1:sz(1),1:sz(2),1:n},sz,sz);

res.Line=obj.Lin;
res.Partition=obj.Par;
res.Slice=obj.Sli;
res.Data=Data;
res.Repetition = obj.Rep;
res.Set = obj.Set;
res.iceParam=obj.iceParam;
res.freeParam=obj.freeParam;
res.timestamp=obj.timestamp;

