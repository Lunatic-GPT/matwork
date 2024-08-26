function recon_DCE(varargin)
% Recon_hahn reads in binary data and displays the corresponding image
% after throwing away the header data.
%
% This .m-file is written for UNIX-OS and not tested for other systems.
%
% Input arguments:
% 1) Image type:
%    - DCE image: varargin{1} = 1
%      --> recon_hahn displays all six images of the DCE series
%    - ADC image: varargin{1} = 2
%      --> recon_hahn displays the ADC image; the same as passing no
%          argument to recon_hahn
%    - DTI image: varargin{1} = 3
%      --> recon_hahn displays all six images of the DTI series
%    - FA image: varargin{1} = 4
%      --> recon_hahn displays the FA image;
% 2) Display properties:
%    - single   slice display: varargin{2} = 0
%    - multiple slice display: varargin{2} = 1
%      --> recon_hahn displays more slices (with respect to chosen slice)
%          > DCE-case: displaying 4 slices inferior and superior
%          > ADC-case: displaying 1 slice inferior and superior
%          > DTI-case: displaying 1 slice inferior and superior
%    - NEW(6/11/08)! NOW DOWNLOAD OF UP TO N SLICES POSSIBLE!: varargin{2} = 2
%      --> Give total number of slices (for DCE) in varargin{8}.
% 3) foldername = varargin{3}
% 4) DTI 'slice-iter', i.e. the 'jump' ('iteration') between the images...
% 5) save_folder
% 6) fnumber
% 7) DCE 'slice-iter'
%
% Created 02/19/2008
% Michigan State University
% Tobias Hahn, on the basis of a reconstruction .m-file provided by Dr.Jie Huang, Michigan State University

%% Parameters
global total_size nf np header_size fname foldername folder fnumber fnum_str save_folder image_type% variables needed for the subfuntion
if nargin < 3
    foldername=input('Enter the folder name of the file of interest: ','s'); % file name of the prime image
    folder = foldername;
else
    folder = varargin{3};
end
k = strfind(folder, 'DCE');
if isempty(k)~=1
    image_type = 'DCE';
end
k = strfind(folder, 'ADC');
if isempty(k)~=1
    image_type = 'ADC';
end
k = strfind(folder, 'DTI');
if isempty(k)~=1
    image_type = 'DTI';
end
k = strfind(folder, 'FA');
if isempty(k)~=1
    image_type = 'FA';
end
foldername = [folder '/MRDC.'];
if nargin < 6
    fnumber = input('Enter the number of your image of interest: ');
else
    fnumber = varargin{6};
end
fnum_str = num2str(fnumber);
fname = numconv(fnumber);

if nargin < 5
    save_folder = '~/matlab/work_images/';
else
    save_folder = varargin{5};
end

if nargin > 6
    slice_iter_DCE = varargin{7};
end

if nargin > 7
    total_slices = varargin{8};
end

%% Input arguments

if nargin == 0
    nf = 512;
    np=nf;
    name = [foldername fname];
    D = dir(name);
    total_size=D.bytes; % get size of file
    header_size=total_size - 2*nf*np;
    image = read_in(name);
    figure; imshow(image, []);
end

if nargin ~= 0
    if varargin{1} == 1 % DCE images
        nf=512;
        np=nf;
        save_subfolder = [save_folder '/' image_type '_temp'];
        if exist(save_subfolder) ~= 0
            rmdir(save_subfolder,'s'); % remove all data in the folder
        end
        mkdir(save_subfolder); % create new empty subfolder
        if varargin{2} == 1
            fnumber = fnumber - 4;
            total_slices = 3;
        elseif varargin{2} == 0
            total_slices = 1;
            % if varargin{2} == 2 --> total_slices is already given in
            % varargin{8}
        end
        fnumber_copy = fnumber;
        if varargin{2} ~= 2
            for slices = 1:total_slices
                fnumber = fnumber_copy + (slices-1)*4; % take only every fourth DCE image
                m = (slices-1)*6; % index of first image of series of specific slice
                slicenumber = num2str(slices);
                for i = 0:5
                    fno = fnumber + slice_iter_DCE*i;
                    filename = numconv(fno);
                    name = [foldername filename];
                    D = dir(name);
                    total_size=D.bytes; % get size of file
                    header_size=total_size - 2*nf*np;
                    image = read_in(name);
                    if (m+i+1) < 10
                        index = num2str(m+i+1);
                        index = ['0' index]; % to ensure the right ordering of the files
                    else
                        index = num2str(m+i+1);
                    end
                    save_filename = [save_subfolder '/' index];
                    save(save_filename, 'image');
                    save_filename = [save_filename '.tiff'];
                    imwrite(image/4000,save_filename,'tif')
                    if i==0
                        quot1 = image;
                    elseif i==1 % save also the WASH-IN image
                        weight_factor = 1;
                        quot2 = image;
                        quot = zeros(size(image,1),size(image,2));
                        quot(quot1~=0) = (quot2(quot1~=0)-quot1(quot1~=0))./(quot1(quot1~=0)*weight_factor);
                        save_filename = [save_subfolder '/wash_in_slice' slicenumber];
                        image = quot;
                        save(save_filename, 'image');
                        save_filename = [save_filename '.tiff'];
                        quot(quot>4)=4;
                        quot = quot/max(max(quot)); % normalization of the wash-in image
                        imwrite(quot,save_filename,'tif');
                    end
                end


                %             % WASH-OUT
                for i=1:6
                    DCE_number = 6*(slices - 1) + i;
                    DCE_number = num2str(DCE_number);
                    if (6*(slices-1)+i)<10
                        DCE_number = ['0' DCE_number];
                    end
                    if nargin < 5
                        file = ['~hahntobi/matlab/work_images/DCE_temp/' DCE_number];
                    else
                        file = [save_folder '/DCE_temp/' DCE_number];
                    end
                    load(file,'image');
                    DCE(:,:,i) = image;
                end
            
            end
        elseif varargin{2} == 2
            for slices = 1:total_slices
                fnumber = fnumber_copy+(slices-1);  % GIVE NOW NUMBER OF FIRST DCE IMAGE IN 'FNUMBER_ARRAY'!!!!
                m = (slices-1)*6; % index of first image of series of specific slice
                slicenumber = num2str(slices);
                for i = 0:5
                    fno = fnumber + slice_iter_DCE*i;
                    filename = numconv(fno);
                    name = [foldername filename];
                    D = dir(name);
                    total_size=D.bytes; % get size of file
                    header_size=total_size - 2*nf*np;
                    image = read_in(name);
                    if (m+i+1) < 10
                        index = num2str(m+i+1);
                        index = ['0' index]; % to ensure the right ordering of the files
                    else
                        index = num2str(m+i+1);
                    end
                    save_filename = [save_subfolder '/' index];
                    save(save_filename, 'image');
                    save_filename = [save_filename '.tiff'];
                    imwrite(image/4000,save_filename,'tif')

                end


                %             % WASH-OUT
                for i=1:6
                    DCE_number = 6*(slices - 1) + i;
                    DCE_number = num2str(DCE_number);
                    if (6*(slices-1)+i)<10
                        DCE_number = ['0' DCE_number];
                    end
                    if nargin < 5
                        file = ['~hahntobi/matlab/work_images/DCE_temp/' DCE_number];
                    else
                        file = [save_folder '/DCE_temp/' DCE_number];
                    end
                    load(file,'image');
                    DCE(:,:,i) = image;
                end
        
        
               
            end
        end
        disp(['DCE images downloaded to ' save_folder '/DCE_temp']);
    end
   
end

close;


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



