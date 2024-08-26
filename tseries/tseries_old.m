function tseries(ts_fname,cc_fname,thr,rep,subbrik_mask,masktype,vml_name,vmr_name)
% tseries(ts_fname,cc_fname,thr,rep,subbrik_mask,masktype[,vml_name,vmr_name])
% ts_fname:  data set of the ordered time series
% cc_fname: data set of the correlation coefficients
% thr:  cc threshold
% rep:  repititions for each stimulus pattern
% subbrik_mask: what subbricks to include in finding 
%               the commonly activated voxels.
%               an array of 1's and 0's with length equal to the number of
%               subbricks
% masktype: whether to include the voxels mapped to more than one visual area
%           rois. 1: yes, 0: no.
% vml_name/vmr_name:  data sets of the left and right visual area masks.
%                   defaults: 'brainmask_sroi2v_lh+orig' and
%                   'brainmask_sroi2v_rh+orig'

if ~exist('vml_name','var')
    vml_name = 'brainmask_sroi2v_lh+orig';
end
if ~exist('vmr_name','var')
    vmr_name = 'brainmask_sroi2v_rh+orig';
end

[maskl,info] = BrikLoad(vml_name);
rois_l = roi_search(info.BRICK_LABS);
[maskr, info] = BrikLoad(vmr_name);
rois_r = roi_search(info.BRICK_LABS);

nrois = length(rois_l);

if (nrois ~=length(rois_r))
    error('The left and rigth hemispheres have different rois');
end

data = BrikLoad(ts_fname);
[cc,info] = BrikLoad(cc_fname);

%% do not modify the code below.
cc_mask = zeros(size(cc));
npatterns = size(cc,4);  % number of patterns
% clusterize the masks
sz = size(cc);
cc_mask_common = ones(sz(1:3));
for i=1:npatterns
cc_mask(:,:,:,i) = clusterize(cc(:,:,:,i)>thr,20,info);
if (subbrik_mask(i))
  cc_mask_common = cc_mask(:,:,:,i) & cc_mask_common;
end

end

umapl = (maskl(:,:,:,1)== 5);  % only map to a single roi.
umapr = (maskr(:,:,:,1)== 5);

figure; 

ind_plot = 1;
subplot(ceil(nrois/2),2,ind_plot);

vmask = zeros([sz(1:3),nrois]);

ts = zeros(size(data,4)/rep,nrois);
er = zeros(size(data,4)/rep,nrois);
for i=1:nrois
if masktype == 0  % only include voxels map to a single roi
    if i<nrois
    vmask(:,:,:,i) = (umapl & maskl(:,:,:,2) ==i) ...
                   | (umapr & maskr(:,:,:,2)==i);
    else
     vmask(:,:,:,i) = maskl(:,:,:,1)== 2 |maskr(:,:,:,1)== 2;    % outside any roi
    end
else  % include voxels map to multiple rois. 
    % Assign such voxels to the roi on which it has the largest number of
    % nodes.
     vmask(:,:,:,i) = (maskl(:,:,:,3+nrois*2+i) == 1 | maskr(:,:,:,3+nrois*2+i) == 1);
end

vi = find(vmask(:,:,:,i) & cc_mask_common );
[ts(:,i), er(:,i)] = add_tseries(data,vi,rep,npatterns);

subplot(ceil((1+nrois)/2),2,i);
errorbar(ts(:,i),er(:,i));
text(20,0.01,rois_l{i},'FontSize',18);

end

brain = ones(size(maskl,1),size(maskl,2),size(maskl,3));
vi = find(brain & cc_mask_common );
[ts_brn,er_brn] = add_tseries(data,vi,rep,npatterns);
subplot(ceil((1+nrois)/2),2,nrois+1);errorbar(ts_brn,er_brn);
text(20,0.01,'brain','FontSize',18);

subbrik_str = sprintf('%d',subbrik_mask);
prefix = strrm(cc_fname,'+orig');

clear data maskl maskr;
fname = sprintf('%s_mask%dcc%3.2fb%s.mat',prefix,masktype,thr,subbrik_str);
save(fname);

function [ts er] = add_tseries(data,vi,rep,ntrl)     

%rep:  reptitions     
     sz = [size(data,1),size(data,2),size(data,3)];
    
     duration = size(data,4)/rep/ntrl;
     ts_tmp = zeros(size(data,4)/rep,rep);     
     
     for i=1:length(vi)
    
      [i1,i2,i3]=ind2sub(sz,vi(i));
      for j=1:rep
          for itr = 1:ntrl
              fst = (itr-1)*duration*rep+(j-1)*duration +1;
              lst = fst + duration -1;
            ts_tmp((itr-1)*duration+1:itr*duration,j) = ...
            ts_tmp((itr-1)*duration+1:itr*duration,j)+ squeeze(data(i1,i2,i3,fst:lst));
          end
      end
     end

     
    ts = mean(ts_tmp,2)/length(vi);
    er = std(ts_tmp,0,2)/sqrt(rep)/length(vi);
      
    
    
function rois=roi_search(label)
         
tilde = find(label=='~');

ntd = length(tilde);

if (mod(ntd-3,3)~=0) 
    error('Wrong roi label format');
end
nrois = (ntd-3)/3;

fprintf('%d rois found\n',nrois);
rois = cell(1,nrois);
for i=1:nrois

    len = 0;
    while label(tilde(3+i)-len-1)~=' ' ...
       && label(tilde(3+i)-len-1)~='/' ...
       && label(tilde(3+i)-len-1)~='\'
      len=len+1;
    end
    
    rois{i}= label(tilde(3+i)-len:tilde(3+i)-1);
    fprintf('%s  ',rois{i});
end



fprintf('\n');



      
         
         
         
         


