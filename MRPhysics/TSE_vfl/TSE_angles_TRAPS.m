function FA_Array=TSE_angles_TRAPS(lKSpaceCenter,lETL)

necho_pss=4;
if lKSpaceCenter>=13
  dpeak_skf = 0.96;
else
  dpeak_skf = 0.92;
end

%lKSpaceCenter = 7;
dpeak=180;
dPI=3.141592654;
protFA=120;  %flip angle from the protocol.



 FA_Array(1) = 90.0 + protFA/2.0 + 0.04*(90.0-protFA/2.0);	%// MW: first refocusser
    PH_Array(1) = 0.0;
    
 %   ///////////////////////////////////////////////////////////////
 %   // ... and continue smooth transition (pseudo steady state) to target flipangle
 %   //     (fill full array, then we are save if something changes in the following implementation)
    drf = FA_Array(1)-protFA;
    
    for e_count=1:lETL-1
    
        FA_Array(e_count+1) = protFA + drf*2.0^(-e_count-0.5);%			// MW: the following refocussers
        
    end


if (lKSpaceCenter > necho_pss)											% MW: 4 x optimized only good for higher KSpaceCenters
                                                                      % MW: between 6 and 8 reduce optimized pulses and go up earlier
            lecho1 = min([necho_pss,lKSpaceCenter-necho_pss]);
            lecho2 = lKSpaceCenter;
            
            drf = dpeak - FA_Array(lecho1);							% MW: NEW_TRAPS sinuosidal up-ramp from FA to peak
            temp_var = lecho2 - lecho1+1;
            lI=1;
            for e_count=lecho1:lecho2-1											
                FA_Array(e_count+1) = dpeak - drf * cos(dPI/2.0 * (lI)/temp_var);
         %       PH_Array(e_count+1) = 0.0;
                lI=lI+1;
            end
else
        % KSCenter <= 4 : no space for ramping -> use flat plateau to KSCenter
            lecho1 = 0;                                                 % MW: no optimized pulses and no ramp, flat beginning is necessary	
            lecho2 = 0;																																
            for e_count=0:lKSpaceCenter-1	        % MW: NEW_TRAPS skaled flat peak
                FA_Array(e_count+1) = dpeak*dpeak_skf;
           %     PH_Array(e_count+1) = 0.0;
      
            end
end
        
       % ///////////////////////////////////////////////////////////////
       % // set peak in KSpaceCenter
       % //
        FA_Array(lKSpaceCenter+1) = dpeak*dpeak_skf;                    %// Note: peak in KSCenter is reduced, 
      %  PH_Array(lKSpaceCenter+1) = 0.0;                                %//       but neihbouring end/start points of ramps are 



      %  ///////////////////////////////////////////////////////////////
       % // ramp down to final plateau
       % //
        lecho3 = lKSpaceCenter+1;
        lecho4 = 2*(lKSpaceCenter+1)-lecho1;                             %   // MW: ramps are basically symmetric


 
        if (lecho4-lecho3 <5)
            lecho4 = lecho3+5;				    %	// MW: ensure minimum down ramp safety   (case: low Isocenter)
        elseif (lecho4-lecho3 >8)
            lecho4 = lecho3+8;					%    // MW: limitation of length of down-ramp (case: long ETL, high Isocenter), good for quality and SAR
        end
        
        drf = dpeak - protFA;
        temp_var = lecho4 - lecho3;                                     %// MW: these two lines HAVE TO BE BEFORE the next one

        if (lecho4 > lETL)
            lecho4 = lETL;
        end %// MW: limit to valid range, BUT do rest of a smooth down ramp with "correct" slope
         lI=1;
        for e_count=lecho3:lecho4-1		%		// MW: NEW_TRAPS sinuosidal down-ramp from peak to end
      							
            FA_Array(e_count+1)= dpeak - drf * sin(dPI/2.0 * lI/temp_var);
            lI=lI+1;
        end

        %//////////////////////////////////////////////////////////////////
        %// flat final flip to end of echotrain
        for e_count=lecho4:lETL-1		       % 	// MW: AUTO_TRAPS flat end
        % // MW: acts only if told (depends on lecho4)
            FA_Array(e_count+1) = protFA;
         
         end
        
        
        
        
        
        
        
        