function writeRFfile_Siemens(rf,fname)

rf = rf./max(abs(rf));

fid = fopen(fname,'w');

fprintf(fid,sprintf('PULSENAME: %s \n',fname));
fprintf(fid,'COMMENT: This is a test file Comment\n');
fprintf(fid,'REFGRAD: 3.70000000 \n');
fprintf(fid,'MINSLICE: 0.200000000 \n');
fprintf(fid,'MAXSLICE: 50.000000000 \n');
fprintf(fid,['AMPINT: ' num2str(sum(rf)) ' \n']);
fprintf(fid,['POWERINT: ' num2str(sum(abs(rf).^2)) '\n']);
fprintf(fid,['ABSINT: ' num2str(sum(abs(rf))) '\n']);

for k = 1:length(rf)
    fprintf(fid,[num2str(abs(rf(k))) ' ' num2str(angle(rf(k)))  '\n']);
end

fclose(fid);