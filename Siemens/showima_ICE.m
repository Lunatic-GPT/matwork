function d=showima_ICE(fname,msize)

a=fopen(fname,'r');

d=fread(a,inf,'uint16');
d=reshape(d,msize);
d=permute(d,[2,1,3,4]); % matrix column corresponds to COL in ICE if the file was saved with WriteToFile function


figure;imshow(d,[]);
fclose(a);
