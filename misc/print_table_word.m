function print_table_word(rlab,clab,m,sep,format)
% print_table_word(rlab,clab,m,sep)
%sep is cell containing symbols to seperate numbers in m; the number of
%symbols should be the same as number of columns in m;
% add ';' to the symbol to divide m into different columns.

if ~isempty(clab)
    fprintf(';%s',clab{:});
end
fprintf('\n');
m=squeeze(m);

for i=1:size(m,1)
    if ~isempty(rlab)
     fprintf('%s;',rlab{i});
    end
    for j=1:size(m,2)
        format_str=[format{j},'%s'];
     fprintf(format_str, m(i,j),sep{j});
    end
    fprintf('\n');
end





