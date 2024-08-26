function dcm2dir_useName(dpattern)
%dcm2analyze(dname,pattern,nslices,nt)


tic;

for i=1:999
 dir_str=dir(sprintf('*%s.%04d*.IMA',dpattern,i));
 if ~isempty(dir_str)
     mkdir(num2str(i));
     movefile(sprintf('*%s.%04d*.IMA',dpattern,i),num2str(i));
 end
end