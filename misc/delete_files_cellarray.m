function delete_files_cellarray(ca)
for i=1:length(ca)
    delete(ca{i});
end