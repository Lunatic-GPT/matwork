function d_pc=epi_flipro_phasecor(d,doflip)

if ~exist('doflip','var')
    doflip=true;
end

if mod(size(d,2),2)==0  %should be 0
 d(:,3,:,:,:)=[];
else
 d(:,1,:,:,:)=[];
end

if doflip
 d(:,1:2:end,:,:,:)=flip(d(:,1:2:end,:,:,:),1);
end

d_pc=epi_phase_corr(d(:,3:end,:,:,:),(d(:,[1,2],:,:,:)));

if mod(size(d_pc,2),2)==0
 d_pc=flip(d_pc,1);   
end