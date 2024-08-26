function B0=B0map_phase(scan_ph,scan_mag)
%B0map(fname,tau,ref)

ph=ri(scan_ph);
ph=ph/(max(ph(:))-min(ph(:)))*2*pi;

mag=ri(scan_mag);


te=readdPar(scan_ph,'EchoTime',true);
te=cell2mat(te);
te=unique(te);

d=mag.*exp(1i*(ph));


for i=1%:length(te)-1
  
    B0(:,:,:,i)=-angle(d(:,:,:,i).*conj(d(:,:,:,i+1)))/(te(i+1)-te(i))/2/pi*1000;
    
end

B0=mean(B0,4);
    
save(sprintf('B0map_%s',scan_ph),'B0');




