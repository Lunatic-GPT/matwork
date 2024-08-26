function DCE_cluster_analysis_for_BLAT(params)

col = get(params, 'left cut-off');
cor = get(params, 'right cut-off');
first_slice = get(params, 'First slice');
nslices = get(params, 'Total number of slices');
DCE_image_dir = get(params, 'Select image directory');
data_dir = get(params, 'Select data directory');
%roi3d = get(params, 'use 3d ROI');  % 
roi3d = false; % disable roi3d for the moment.

dir_struct = dir(fullfile(data_dir,'*X*'));
if isempty(dir_struct)
    error(['Wrong data directory: ', data_dir]);
end

set(gcbo,'Enable','off');
results = {};
%% load DCE images
h_w = waitbar(0,'Detailed analysis');

DCE = zeros(512,512,nslices,6);

for slice = 1:nslices
    for i=1:6 % six DCE images
        file = fullfile(DCE_image_dir,sprintf('%02d_%d_shifted.mat',first_slice+slice-1,i));
        load(file,'image');
        DCE(:,:,slice,i) = double(image);
    end
end

%% Load wash in
waitbar(1/8,h_w);
wash_in_nd = zeros(512,512,nslices);
for slice = 1:nslices
    file = fullfile(DCE_image_dir, sprintf('%02d_wash_in_norm_diff.mat',slice+first_slice-1) ); 
    load(file,'image');
    wash_in_nd(:,:,slice) = image;
end

%% Load rois
waitbar(2/8,h_w);
DCE_inner = zeros(512,512,nslices);
DCE_tissue = zeros(512,512,nslices);
DCE_tissue2 = zeros(512,512,nslices);

if ~roi3d
  for slice = 1:nslices
    dir_struct = dir(fullfile(data_dir,['*X',num2str(slice+first_slice-1)]));
    
    if(length(dir_struct)>1)
        disp(['Multiple directories for slice ',num2str(slice), ' in ',data_dir]);
        disp(['Use ',dir_struct(1).name]);
        dir_struct = dir_struct(1);
    end
    
    [token,rem] = strtok(dir_struct.name, 'X');
    if isempty(rem)
        token = rem;
    end
    DCE_dir = fullfile(data_dir,[token 'X' num2str(slice+first_slice-1)]);
    
    final_BW_file = fullfile(DCE_dir, 'ROIs.mat');
    load(final_BW_file, 'inner_BW_large', 'tissue_BW_large', 'tissue2_BW_large');
    DCE_inner(:,:,slice) = inner_BW_large;
    DCE_tissue(:,:,slice) = tissue_BW_large;
    DCE_tissue2(:,:,slice) = tissue2_BW_large;
  end
else
    final_BW_file = fullfile(data_dir,'3dROI', 'ROI3d.mat');
    load(final_BW_file, 'inner_BW_large', 'tissue_BW_large', 'tissue2_BW_large');
    DCE_inner = inner_BW_large(:,:,first_slice:first_slice+nslices-1);
    DCE_tissue = tissue_BW_large(:,:,first_slice:first_slice+nslices-1);
    DCE_tissue2 = tissue2_BW_large(:,:,first_slice:first_slice+nslices-1);
end
    
    
    
%%  signal analysis in roi, tissue, and tissue 2.
waitbar(3/8,h_w);
DCEimage_tissue_vec = zeros(512,512,nslices,6);
DCEimage_tissue2_vec = zeros(512,512,nslices,6);
DCEimage_lesion_vec = zeros(512,512,nslices,6);

 
for i=1:6
    
  DCEimage_tissue_vec(:,:,:,i) = DCE(:,:,:,i).*DCE_tissue;
  DCEimage_tissue2_vec(:,:,:,i) = DCE(:,:,:,i).*DCE_tissue2;
  DCEimage_lesion_vec(:,:,:,i) = DCE(:,:,:,i).*DCE_inner;
end

 
 for i=1:6
     temp = DCEimage_tissue_vec(:,:,:,i);
     temp = temp(:);
     DCEimage_tissue_mean(i) = mean(temp(temp~=0));
     DCEimage_tissue_std(i) = std(temp(temp~=0));
 
     temp = DCEimage_tissue2_vec(:,:,:,i);
     temp = temp(:);
     DCEimage_tissue2_mean(i) = mean(temp(temp~=0));
     DCEimage_tissue2_std(i) = std(temp(temp~=0));
 
     temp = DCEimage_lesion_vec(:,:,:,i);
     temp = temp(:);
     DCEimage_lesion_mean(i) = mean(temp(temp~=0));
     DCEimage_lesion_std(i) = std(temp(temp~=0));
 end   
            
results(end+1,1:2) = {DCEimage_lesion_mean 'DCE image no.1-6 lesion mean'}; %
results(end+1,1:2) = {DCEimage_tissue_mean 'DCE image no.1-6 tissue mean'}; %
results(end+1,1:2) = {DCEimage_tissue2_mean 'DCE image no.1-6 tissue2 mean'}; %

results(end+1,1:2) =  {DCEimage_lesion_std 'DCE image no.1-6 lesion std'}; %
results(end+1,1:2) = {DCEimage_tissue_std 'DCE image no.1-6 tissue std'}; %
results(end+1,1:2) = {DCEimage_tissue2_std 'DCE image no.1-6 tissue2 std'}; %

%% WASH-IN analysis with DCE ROI
waitbar(4/8,h_w);
   wash_in_nd_ROI = wash_in_nd.*DCE_inner;
   temp = wash_in_nd_ROI(:);
   wash_in_nd_ROI_DCE_mean = mean(temp(temp~=0));
   wash_in_nd_ROI_DCE_std = std(temp(temp~=0));
    

    wash_in_nd_tissue = wash_in_nd.*DCE_tissue;
    temp = wash_in_nd_tissue(:);
    wash_in_nd_tissue_mean = mean(temp(temp~=0));
    wash_in_nd_tissue_std = std(temp(temp~=0));
    
    wash_in_nd_tissue2 = wash_in_nd.*DCE_tissue2;
    temp = wash_in_nd_tissue2(:);
    wash_in_nd_tissue2_mean = mean(temp(temp~=0));
    wash_in_nd_tissue2_std = std(temp(temp~=0));
    
results(end+1,1:2) = {wash_in_nd_ROI_DCE_mean, 'wash-in norm diff lesion mean'}; %
results(end+1,1:2) = {wash_in_nd_ROI_DCE_std,'wash-in norm diff lesion std' }; %
results(end+1,1:2) = {wash_in_nd_tissue_mean, 'wash-in norm diff tissue mean'}; %
results(end+1,1:2) =  {wash_in_nd_tissue_std 'wash-in norm diff tissue std' }; %
results(end+1,1:2) =  {wash_in_nd_tissue2_mean 'wash-in norm diff tissue 2 mean'}; %
results(end+1,1:2) = {wash_in_nd_tissue2_std 'wash-in norm diff tissue 2 std'}; %


%% WASH-OUT Analysis
waitbar(5/8,h_w);
wash_out = zeros(512,512,nslices);
for slice = 1:nslices
    file = fullfile(DCE_image_dir, sprintf('%02d_wash_out.mat', slice+first_slice-1));
    
    load(file,'image');
    wash_out(:,:,slice) = image;
end

    wash_out_ROI = wash_out.*DCE_inner;
    
    temp = wash_out_ROI(:);
    wash_out_ROI_DCE_mean = mean(temp(temp~=0));
    wash_out_ROI_DCE_std = std(temp(temp~=0));


    wash_out_tissue = wash_out.*DCE_tissue;
    temp = wash_out_tissue(:);
    wash_out_tissue_mean = mean(temp(temp~=0));
    wash_out_tissue_std = std(temp(temp~=0));
  
    wash_out_tissue2 = wash_out.*DCE_tissue2;
    temp =  wash_out_tissue2(:);
    wash_out_tissue2_mean = mean(temp(temp~=0));
    wash_out_tissue2_std = std(temp(temp~=0));
    
results(end+1,1:2) =  {wash_out_ROI_DCE_mean 'wash-out lesion mean'}; %
results(end+1,1:2) = {wash_out_ROI_DCE_std 'wash-out lesion std'}; %
results(end+1,1:2) = {wash_out_tissue_mean 'wash-out tissue mean'}; %
results(end+1,1:2) = {wash_out_tissue_std 'wash-out tissue std'}; %
results(end+1,1:2) = {wash_out_tissue2_mean 'wash-out tissue 2 mean'}; %
results(end+1,1:2) = {wash_out_tissue2_std 'wash-out tissue 2 std'}; %
%% lesion volume


 results(end+1,1:2) = {sum(DCE_inner(:)), 'lesion volume'};
 
 
%% NEW WASH-IN (in degrees): created 6/28/2008
waitbar(6/8,h_w);
sz = size(DCE); 
wash_in = zeros(sz(1:3));

for i=1:nslices
    fname = fullfile(DCE_image_dir,sprintf('%02d_wash_in.mat',i+first_slice-1));
    load(fname,'image');
    wash_in(:,:,i) = image;
    
end 
    wi_lesion = wash_in.*DCE_inner;
    wi_tissue = wash_in.*DCE_tissue;
    wi_tissue2 = wash_in.*DCE_tissue2;
  
   temp = wi_lesion(:);
   wi_lesion_mean = mean(temp(temp~=0));
   temp = wi_tissue(:);
   wi_tissue_mean = mean(temp(temp~=0));
   temp = wi_tissue2(:);
   wi_tissue2_mean = mean(temp(temp~=0));

   slope_angle = (180-wi_lesion+wash_out_ROI).*DCE_inner;
   temp = slope_angle(:);
   
   slope_mean = mean(temp(temp~=0));
    
results(end+1,1:2) = {wi_lesion_mean 'new wi_lesion_mean'}; %
results(end+1,1:2) = {wi_tissue_mean 'new wi_tissue_mean'}; %
results(end+1,1:2) = {wi_tissue2_mean 'new wi_tissue2_mean'}; %
results(end+1,1:2) = {slope_mean 'slope angle'}; %


savename = fullfile(data_dir, 'analysis_results.mat');
save(savename, 'results','params');
disp(['Analysis results saved to ',savename]);

waitbar(7/8,h_w);
results_detailed_analysis = DCE_detailed_ROI_analysis(DCE,wash_in,wash_out,DCE_inner,col,cor);

savename = fullfile(data_dir, 'detailed_analysis_results.mat');
save(savename, 'results_detailed_analysis','params');
disp(['Detailed analysis results saved to ',savename]);
delete(h_w);
disp('Finished');
set(gcbo,'Enable','on');

