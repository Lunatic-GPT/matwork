function [Trange,Rrange,MT,MR]=motionRange_dicom(d_moco,i_start,i_end)

tmp=Moco_ImageComments(d_moco);

if ~exist('i_start','var')
    i_start=1;
end

if ~exist('i_end','var')
    i_end=size(tmp,1);
end


pmin=min(tmp(i_start:i_end,:),[],1);
pmax=max(tmp(i_start:i_end,:),[],1);

Trange=pmax(4:6)-pmin(4:6);
Rrange=pmax(1:3)-pmin(1:3);

MR=sos(pmin(1:3)-pmax(1:3));
MT=sos(pmin(4:6)-pmax(4:6));



