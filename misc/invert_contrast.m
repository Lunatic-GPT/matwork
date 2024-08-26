function invert_contrast(fpattern)

list=dir(fpattern);

for i=1:length(list)
    
   disp(list(i).name);
   prefix=strtok(list(i).name,'.');

   fname=fullfile(list(i).folder,list(i).name);

   save_name=filename_append(list(i).name,'_Inv',1);
    
   if exist(save_name,'file')
       continue;
   end

    d=load_untouch_niigz(fname);
    
    d.img=max(d.img(:))-d.img;
    save_untouch_niigz(d,save_name);

    
end
