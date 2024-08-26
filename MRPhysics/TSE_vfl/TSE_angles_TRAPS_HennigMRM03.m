function FA_Array=TSE_angles_TRAPS_HennigMRM03(n,lETL,a1,a4)


  dpeak_skf = 0.92;

%lKSpaceCenter = 7;
dpeak=180;
dPI=3.141592654;
protFA=a1;  %flip angle from the protocol.

 FA_Array(1) = 90.0 + protFA/2.0 + 0.04*(90.0-protFA/2.0);	%// MW: first refocusser
  
 %   ///////////////////////////////////////////////////////////////
 %   // ... and continue smooth transition (pseudo steady state) to target flipangle
 %   //     (fill full array, then we are save if something changes in the following implementation)
    drf = FA_Array(1)-protFA;
    
    for e_count=1:n(1)-1
    
        FA_Array(e_count+1) = protFA + drf*3.5^(-e_count-0.5);%			// MW: the following refocussers
        
    end

    tmp=linspace(a1,180,n(2)-n(1)+1);
FA_Array(end+1:end+n(2)-n(1))=tmp(2:end);

if n(3)>n(2)
  FA_Array(end+1:end+n(3)-n(2))=180;
end

    tmp=linspace(180,a4,n(4)-n(3)+1);
FA_Array(end+1:end+n(4)-n(3))=tmp(2:end);

if lETL>n(4)
  FA_Array(end+1:end+lETL-n(4))=a4;
end

    
figure;plot(FA_Array);
        
        
        
        
        