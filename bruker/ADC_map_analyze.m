function ADC_map_analyze(prefix)
  [a,d]=readanalyze(prefix);
  
b=readbPar([prefix,'/method'],'PVM_DwEffBval');

adc=log(mean(a(:,:,:,b<50),4)./mean(a(:,:,:,b>50),4));

adc=adc/mean(b(b>50))*1000;
writeanalyze(adc,[prefix,'_adc'],d);



