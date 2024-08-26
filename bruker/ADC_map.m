function ADC_map(prefix,irecon)
if ~exist('irecon','var')||irecon==1
  a=BrikLoad([prefix,'+orig']);
else
    a=BrikLoad(sprintf('%s_%d+orig',prefix,irecon));
end

b=readbPar([prefix,'/method'],'PVM_DwEffBval');

adc=log(a(:,:,:,b<50)./mean(a(:,:,:,b>50),4));

adc=adc/mean(b(b>50))*1000;
write_afni(adc,[prefix,'_adc']);



