function recon_DCE(source_folder,save_folder,fslice,DCE_iter,nslices)
%% Recon_hahn reads in binary data and displays the corresponding image
% after throwing away the header data.
%
% This .m-file is written for UNIX-OS and not tested for other systems.
%
% Input arguments:
% 1) source_folder:   the source folder containing the dicom images
% 2) save_folder: the folder to save the images
% 4) fnumber: the first slice for analysis. 1 based
% 5) DCE_iter: the total number of slices in each volume
% 6) nslices: total number of slices to analyze.
%
% Created 02/19/2008
% Michigan State University
% Tobias Hahn, on the basis of a reconstruction .m-file provided by Dr.Jie Huang, Michigan State University
% 4-28-2009:  remove DTI part, XZ.

%% Parameters
global total_size nf np header_size;% variables needed for the subfuntion



    nf = 512;
    np=nf;
    

    
    
        if exist(save_folder,'dir')
            rmdir(save_folder,'s'); % remove all data in the folder
        end

        mkdir(save_folder); % create new empty subfolder
    %    h=waitbar(0,'download images...','Units','normalized','Position',[0.3,0.5,0.12,0.05]);
        h=waitbar(0,'download images...');
           for slices = 1:nslices
                fnumber = fslice+(slices-1);  %
               
                         
                for i = 0:5
                    fno = fnumber + DCE_iter*i;
                    
                    name = sprintf('%s/MRDC.%04d',source_folder, fno);
                    D = dir(name);
                    total_size=D.bytes; % get size of file
                    header_size=total_size - 2*nf*np;
                    image = read_in(name);
                  
                    index = sprintf('%02d_%d',slices,i+1);
                    save_filename = fullfile(save_folder, index);
                    save(save_filename, 'image');
                    imwrite(image/4000,[save_filename,'.tiff'],'tif')
               
                end      
                waitbar(slices/nslices,h);
            end
             
            
    disp(['DCE images downloaded to ' save_folder]);

  delete(h);  
%close;


%%
function image = read_in(filename)
global nf np header_size 
% Read in of file
image=zeros(np,nf);
offset=header_size; % given in bytes
fid=fopen(filename,'r');
fseek(fid,offset,'bof'); % repositions the pointer right after the header data
raw_size=nf*np;
[rawdat,COUNT]=fread(fid,raw_size,'int16'); 
fclose(fid);
if COUNT ~= raw_size
 warning('There is an unknown error with the image size');
end

for j=1:np %phase  
  for i=1:nf %frequency
    image(j,i)=rawdat(i+(j-1)*nf);
  end
end



