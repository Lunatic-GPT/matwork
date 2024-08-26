function FOVcenter(a2)


oa=readbPar([a2,'/method'],'PVM_SPackArrSliceOffset');


na=readbPar([a2,'/acqp'],'NSLICES');
dist=readbPar([a2,'/method'],'PVM_SPackArrSliceDistance');

fprintf('offset = %f; nslice = %d; thickness = %f\n',oa,na,dist);
