function [adc,adc2]=tgse_vfl_pss_ADCEnable(etl,epi,stb,lLinesToMeasure,lPartitionsToMeasure)
% adc=tgse_vfl_pss_ADCEnable(etl,epi,stb,lLinesToMeasure,lPartitionsToMeasure)
% etl: ETL, do not include the effect of stb.
% epi: epi factor for each spin echo.
% stb: slice turbo factor

lRepetitions=floor((lPartitionsToMeasure+stb-1)/stb);

adc=zeros(etl*epi*stb,lRepetitions);

lin=adc;
par=adc;

adc2=zeros(etl*epi*stb,lRepetitions);


for lRep=0:lRepetitions-1
    for lI=0:etl-1
        for lJ=0:stb-1
            
            offset=lJ*epi+lI*epi*stb+lRep*epi*stb*etl;
            
            for lK=0:epi-1
                lLinNo=lI+lK*etl;
                
                if (lLinNo<lLinesToMeasure)
                    adc(lK+offset+1)=1;
                    lin(lK+offset+1)=lLinNo;
                else
               %     adc(lK+offset+1)=0;
                    lin(lK+offset+1)=lLinesToMeasure-1;
                end
                
                lParNo = lJ+lRep*stb;
                
                
                if (lParNo<lPartitionsToMeasure)
            %        adc(lK+offset+1)=0;
                    par(lK+offset+1)=lParNo;
                else
                    adc(lK+offset+1)=0;
                    par(lK+offset+1)=lPartitionsToMeasure-1;
                end
            end
        end
    end
end



for lRep=0:lRepetitions-1
    for lI=0:etl-1
        for lJ=0:stb-1
            
            offset=lJ*epi+lI*epi*stb+lRep*epi*stb*etl;
            
            offset2=lRep*epi*stb*etl;
            
            for lK=0:epi-1
                lADCLoopCounter = lK+(lJ+lI*stb)*epi;                
                lEchoLoopCounter=(lJ+lI*stb);
                res = lADCLoopCounter<epi*stb;
                
                if ~res
                    test2 = lin(lADCLoopCounter+1+offset2)~=lin(lADCLoopCounter+1+offset2-epi*stb);
                    test3 = par(lADCLoopCounter+1+offset2)~=par(lADCLoopCounter+1+offset2-epi);
                    test4=mod(lEchoLoopCounter,stb)==0;
                    res=test2&&(test4||test3);
                end
                
                adc2(lK+offset+1)=res;
                
            end
        end
    end
end



