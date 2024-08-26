function zeroslice_FOVmatch(a,a2,b)

[data,dim]=readanalyze(a);
da=readbPar([a2,'/method'],'PVM_SPackArrSliceDistance');

db=readbPar([b,'/method'],'PVM_SPackArrSliceDistance');


if da~=db
    error('da ~= db');
end
oa=readbPar([a2,'/method'],'PVM_SPackArrSliceOffset');
ob=readbPar([b,'/method'],'PVM_SPackArrSliceOffset');


na=readbPar([a2,'/acqp'],'NSLICES');
nb=readbPar([b,'/acqp'],'NSLICES');


posa_min=oa-(na-1)/2*da;
posb_min=ob-(nb-1)/2*da;

posa_max=oa+(na-1)/2*da;
posb_max=ob+(nb-1)/2*da;

sz=size(data);
ns=(posa_min-posb_min)/da;
    
if ns>0    
    data=cat(3,zeros([sz(1:2),ns]),data);
else
    data(:,:,1:abs(ns))=[];
end


ns=(posb_max-posa_max)/da; 
if ns>0    
    data=cat(3,data,zeros([sz(1:2),ns]));
else
    data(:,:,end-abs(ns)+1:end)=[];
end

a=strtok(a,'.');

writeanalyze(data,[a,'_zf'],dim);

