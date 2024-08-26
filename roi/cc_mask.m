function cc_mask(ccfile,thr,pats)
%cc_mask(ccfile,thr,pats)
brik_sel = [];
expr = [];
prefix = [];
for i=1:length(pats)
    
brik_sel = sprintf('%s -%c %s''[%d]''',brik_sel,'a'+i-1,ccfile,2*pats(i)-2);   

prefix = sprintf('%s%d_',prefix,pats(i));
if i<length(pats)
expr = sprintf('%sstep(%c-%3.2f)*',expr,'a'+i-1,thr);
else
    expr = sprintf('%sstep(%c-%3.2f)',expr,'a'+i-1,thr);
end
end

ccfile_pre = strtok(ccfile,'+');
cmd = sprintf('3dcalc -prefix %s_pat_%sthr%3.2f %s -expr ''%s''', ...
               ccfile_pre,prefix,thr,brik_sel,expr);
           unix(cmd);