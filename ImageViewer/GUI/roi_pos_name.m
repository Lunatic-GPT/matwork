function name = roi_pos_name(slice,id,isc)
% allow negative slices
% isc: if true the name will contain c at the end

if ~exist('isc','var')
    isc=0;
end

if isc==0
name=sprintf('roi_pos%d_%d',slice,id);
else
    name=sprintf('roi_pos%d_%d_%d',slice,id,isc);
end
name=strrep(name,'-','_');


    

   