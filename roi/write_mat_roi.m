function m=write_mat_roi(m,prefix)

       
       fo
           
for i=1:size(m,3)
    if any(vec(res(:,:,i))>0)
        roi=m(:,:,i);
        save(sprintf('%sX%d.mat',prefix,i),'roi');
    end
end

       
       
       
