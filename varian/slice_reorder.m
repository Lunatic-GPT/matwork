function data=slice_reorder(data,pss)
%slice_reorder(data,pss)
   if pss ==1
       return;
   end
   
   [pss,ind] = sort(pss);
   
   data=data(:,:,ind,:);
