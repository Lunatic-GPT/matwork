function writeRFfile(rf)

rf = rf./max(rf);

fid = fopen('SLR.90_2mm.pta','w');

fprintf(fid,'PULSENAME: SLR.90_2mm.pta \n');
fprintf(fid,'COMMENT: This is a test file Comment\n');
fprintf(fid,'REFGRAD: 3.70000000 \n');
fprintf(fid,'MINSLICE: 2.000000000 \n');
fprintf(fid,'MAXSLICE: 2.000000000 \n');
fprintf(fid,['AMPINT: ' num2str(sum(rf)) ' \n']);
fprintf(fid,['POWERINT: ' num2str(sum(abs(rf).^2)) '\n']);
fprintf(fid,['ABSINT: ' num2str(sum(abs(rf))) '\n']);

for k = 1:length(rf)
    fprintf(fid,[num2str(abs(rf(k))) ' ' num2str(angle(rf(k))) '; Sampling point: ' num2str(k) '\n']);
end

fclose(fid);