function tiff_rm_alpha(fname)
%% for tiff images generated in powerpoint, open and save the images in
%% painter before run this program to get correct results.

flist=dir(fname);

for i=1:length(flist)
    fname=flist(i).name;
prefix=strtok(fname,'.');
m=imread([prefix,'.tif']);
if size(m,3)==3
    disp('No alpha channel found');
    return;
end
movefile([prefix,'.tif'],[prefix,'_orig.tif'],'f'); 

mm=repmat(m(:,:,4)==0,[1,1,3]);
m=m(:,:,1:3);

m(mm)=255;

imwrite(m(:,:,1:3),[prefix,'.tif'],'TIFF','Resolution',[300,300]);

end
