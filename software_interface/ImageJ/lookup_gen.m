fid=fopen('jet_whiteTop.lut','w');
%cm=uint16(gray(256)*255);
%cm(256,:)=[0,255,0];



cm=uint16([0,0,0;jet(254)*254;255,255,255]);

fwrite(fid,cm(:,1));
fwrite(fid,cm(:,2));
fwrite(fid,cm(:,3));

fclose(fid);


%%

