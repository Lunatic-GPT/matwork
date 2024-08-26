function results = DCE_detailed_ROI_analysis(DCE,wash_in,wash_out,DCE_inner,col,cor)



results = {}; % result - cell


      wash_out = wash_out.*DCE_inner;
    
      BW_wo = (wash_out<col & wash_out>-90).*DCE_inner;
      BW_pers = (wash_out>cor & wash_out<90).*DCE_inner;
      BW_plat = (wash_out>col & wash_out <cor).*DCE_inner;
   
%% image intensity 
    
    
    for i=1:6
        temp = DCE(:,:,:,i).*BW_pers;
        temp = temp(:);
        DCEimage_pers_mean(i) = mean(temp(temp~=0));
        DCEimage_pers_std(i) = std(temp(temp~=0));
        
        temp = DCE(:,:,:,i).*BW_plat;
        temp = temp(:);
        DCEimage_plat_mean(i) = mean(temp(temp~=0));
        DCEimage_plat_std(i) = std(temp(temp~=0));
   
        temp = DCE(:,:,:,i).*BW_wo;
        temp = temp(:);
        DCEimage_wo_mean(i) = mean(temp(temp~=0));
        DCEimage_wo_std(i) = std(temp(temp~=0));
       
    end
    
    results(end+1,1:2) = {DCEimage_wo_mean 'DCE image no.1-6 wo mean'};
    results(end+1,1:2) = {DCEimage_pers_mean 'DCE image no.1-6 pers mean'};
    results(end+1,1:2) = {DCEimage_plat_mean 'DCE image no.1-6 plat mean'}; 

    results(end+1,1:2) = {DCEimage_wo_std 'DCE image no.1-6 wo std'}; 
    results(end+1,1:2) = {DCEimage_pers_std 'DCE image no.1-6 pers std'};
    results(end+1,1:2) = {DCEimage_plat_std 'DCE image no.1-6 plat std'};

%% relative changes
  DCE_ref = DCE(:,:,:,1);
  DCE_ref(DCE_ref==0)=1; % setting zeros to one, we will forget these pixels later by using the 'no_zero_mask', just set so that the division by DCE_ref is possible
    
    for i=2:6
        
        relative_change = (DCE(:,:,:,i)-DCE_ref)./DCE_ref;
        
        temp = relative_change.*BW_wo;
        temp = temp(:);
        relative_change_wo_mean(i) = mean(temp(temp~=0));
        relative_change_wo_std(i) = std(temp(temp~=0));
        

        temp = relative_change.*BW_pers;
        temp = temp(:);
        relative_change_pers_mean(i) = mean(temp(temp~=0));
        relative_change_pers_std(i) = std(temp(temp~=0));
        

        temp = relative_change.*BW_plat;
        temp = temp(:);
        relative_change_plat_mean(i) = mean(temp(temp~=0));
        relative_change_plat_std(i) = std(temp(temp~=0));
        
    end
        
    results(end+1,1:2) =  {relative_change_wo_mean 'relative change wo image i-1 mean'}; 
    results(end+1,1:2) = {relative_change_wo_std 'relative change wo image i-1 std'}; 
    results(end+1,1:2) = {relative_change_plat_mean 'relative change plat image i-1 mean'};
    results(end+1,1:2) = {relative_change_plat_std 'relative change plat image i-1 std'};
    results(end+1,1:2) = {relative_change_pers_mean 'relative change pers image i-1 mean'};
    results(end+1,1:2) = {relative_change_pers_std 'relative change pers image i-1 std'};

%% WASH-IN Analysis
    temp = wash_in.*BW_wo;
    temp = temp(:);
    wash_in_wo_mean = mean(temp(temp~=0));
    wash_in_wo_std = std(temp(temp~=0));
        
    temp = wash_in.*BW_pers;
    temp = temp(:);
    wash_in_pers_mean = mean(temp(temp~=0));
    wash_in_pers_std = std(temp(temp~=0));
        
    temp = wash_in.*BW_plat;
    temp = temp(:);
    wash_in_plat_mean = mean(temp(temp~=0));
    wash_in_plat_std = std(temp(temp~=0));
        
    results(end+1,1:2) = {wash_in_wo_mean 'wash-in wo mean'};
    results(end+1,1:2) = {wash_in_wo_std 'wash-in wo std'};
    results(end+1,1:2) = {wash_in_pers_mean 'wash-in pers mean'};
    results(end+1,1:2) = {wash_in_pers_std 'wash-in pers std'};
    results(end+1,1:2) = {wash_in_plat_mean 'wash-in plat mean'};
    results(end+1,1:2)= {wash_in_plat_std 'wash-in plat std'};

%%    %% WASH-OUT Analysis

    
    temp = wash_out.*BW_wo;
    temp = temp(:); 
    wash_out_wo_mean = mean(temp(temp~=0));
    wash_out_wo_std = std(temp(temp~=0));
   

    temp= wash_out.*BW_pers;
    temp = temp(:);
    wash_out_pers_mean = mean(temp(temp~=0));
    wash_out_pers_std = std(temp(temp~=0));
    
    temp = wash_out.*BW_plat;
    temp = temp(:);
    wash_out_plat_mean = mean(temp(temp~=0));
    wash_out_plat_std = std(temp(temp~=0));
    
results(end+1,1:2) = {wash_out_wo_mean 'wash-out wo mean'};
results(end+1,1:2) = {wash_out_wo_std 'wash-out wo std'};
results(end+1,1:2) = {wash_out_pers_mean 'wash-out pers mean'};
results(end+1,1:2) = {wash_out_pers_std 'wash-out pers std'};
results(end+1,1:2) = {wash_out_plat_mean 'wash-out plat mean'};
results(end+1,1:2) = {wash_out_plat_std 'wash-out plat std'};

%% Volume
    BW_wo_volume = sum(BW_wo(:))/sum(DCE_inner(:));
    BW_plat_volume = sum(BW_plat(:))/sum(DCE_inner(:));
    BW_pers_volume = sum(BW_pers(:))/sum(DCE_inner(:));
 
        results(end+1,1:2) = {BW_wo_volume 'BW_wo_volume fraction'};
        results(end+1,1:2) = {BW_plat_volume 'BW_plat_volume fraction'};
        results(end+1,1:2) = {BW_pers_volume 'BW_pers_volume fraction'};
        
        



