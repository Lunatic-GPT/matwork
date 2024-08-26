function dCBVmap(fname,fname_dtr,bl,te,dr_mion,suffix,dr_bold)
%dCBVmap(fname,fname_dtr,bl,te,dr_mion[,prefix,dr_bold])
% calculate relative CBV changes during stimulation.
% fname: trial mean time course.
% fname_dtr: detrended trial mean time course.
% bl: baseline time points, 1 based.
% te: te during CBV response measurement.
% dr_mion: R2* changes due to MION injection
% prefix: output prefix for normalized dCBV and dR2. Default ''.
% dr_bold: R2* changes due to BOLD.
% dr_bold can be shorter than the CBV time course. However, the stimulus
% presentation time should match.

[b,info]=BrikLoad(fname);
bdtr=BrikLoad(fname_dtr);

s0 = mean(b(:,:,:,bl),4);

b_dtr_n = bdtr./repmat(s0,[1,1,1,size(b,4)]);

dr = -log(b_dtr_n+1)/te;
if ~exist('prefix','var')
    prefix = '';
end
if exist('dr_bold','var')
  drb = BrikLoad(dr_bold);
  l = size(drb,4);
  dr_bc=dr;
  dr_bc(:,:,:,1:l) = dr(:,:,:,1:l)-drb;
  
else
    dr_bold = '';
end


dr0 = BrikLoad(dr_mion);

dCBV = dr./repmat(dr0,[1,1,1,size(dr,4)]);
dCBV_bc = dr_bc./repmat(dr0,[1,1,1,size(dr,4)]);


history = sprintf('dCBVmap(%s,%s,%s,%f,%s,%s,%s)',fname,fname_dtr,num2str(bl),te,dr_mion,prefix,dr_bold);


WriteBrikEZ(dCBV,info,history,['dCBV',suffix]);
WriteBrikEZ(dr,info,history,['dR2',suffix]);

if ~isempty(dr_bold)
  WriteBrikEZ(dCBV_bc,info,history,['dCBV_BC',suffix]);
  WriteBrikEZ(dr_bc,info,history,['dR2_BC',suffix]);
end

