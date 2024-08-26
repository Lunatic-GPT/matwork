function writeRF_epic(m,ph,prefix)
% writeRF_epic(m,ph,prefix)
% m: magnitude
% ph: phase
fid = fopen([prefix,'.rho'],'w','ieee-be');

m = round(m*(2^15-2)/max(m));
fwrite(fid,m,'int16');
fclose(fid);

fid = fopen([prefix,'.pha'],'w','ieee-be');
ph = mod(ph,2*pi)-pi;
ph = round(ph*(2^15-2)/pi);

fwrite(fid,ph,'int16');
fclose(fid);
