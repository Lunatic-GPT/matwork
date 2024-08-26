function tseries_crssvalid(ts_fname,cc_fname,out_fname,thr,rep,npat,lr,subbrik_mask,masktype,vml_name,vmr_name)
% tseries(ts_fname,cc_fname,out_fname,thr,rep,npat,[lr,subbrik_mask,masktype,vml_name,vmr_name])
% ts_fname:  data set of the ordered time series
% cc_fname: data set of the correlation coefficients
% thr:  cc threshold
% rep:  repititions for each stimulus pattern
% subbrik_mask: what subbricks to include in finding 
%               the commonly activated voxels.
%               an array of 1's and 0's with length equal to the number of
%               subbricks.  Default: search all subbricks.
% masktype: whether to include the voxels mapped to more than one visual area
%           rois. 1: yes, 0: no. Default: 1
% vml_name/vmr_name:  data sets of the left and right visual area masks.
%                   defaults: 'brainmask_sroi2v_lh+orig' and
%                   'brainmask_sroi2v_rh+orig'
% lr: lh, rh, or both.
% 4-10-2009:  added lr to allow a single hemisphere.
% 4-27-2009:  added cc_xflt to cc_fname list such that the correlation coefficient 
%             can be calculated from the concatenated time course from different filters  
% 4-30-2009:  allow cc_xflt to be included or not included.

% variables dimensions:  ts: [len_per_trial,rep,nrois,npat,nflt]
%                        cc_mask: [:,:,:,npat,nflt]
%                        cc_mask_xpat: [:,:,:,nflt]
%                        cc_mask_xflt: [:,:,:,npat]
%                        cc_mask_xpatxflt: [n1,n2,n3]
%                        cc: [:,:,:,npat,nflt]
%                        data: [:,:,:,:,nflt]
%
% added mask across both filters and patterns mask_xpatxflt

if ~exist('vml_name','var')
    vml_name = 'brainmask_sroi2v_lh+orig';
end

if ~exist('vmr_name','var')
    vmr_name = 'brainmask_sroi2v_rh+orig';
end

if ~exist('masktype','var')
    masktype = 1;
end

if ~exist('subbrik_mask','var')
    subbrik_mask = ones(1,npat);
end

if ~exist('lr','var')
    lr = 'both';
end

nflt = length(ts_fname);
ncc = length(cc_fname);
if nflt+1 ~= length(cc_fname) && nflt ~= ncc
    error('Incorrect number of corr coef files');
end

[maskl,info] = BrikLoad(vml_name);
rois_l = roi_search(info.BRICK_LABS);
[maskr, info] = BrikLoad(vmr_name);
rois_r = roi_search(info.BRICK_LABS);

nrois = length(rois_l);

if (nrois ~=length(rois_r))
    error('The left and rigth hemispheres have different rois');
end

n1 = size(maskl,1);
n2 = size(maskl,2);
n3 = size(maskl,3);

for iflt =1:nflt 
    
   data(:,:,:,:,iflt) = BrikLoad(ts_fname{iflt});
   [cc(:,:,:,:,iflt),info] = BrikLoad(cc_fname{iflt});

end



cc_mask = zeros(n1,n2,n3,npat,nflt);
cc_mask_xpat = ones(n1,n2,n3,nflt);
cc_mask_xflt = ones(n1,n2,n3,npat);  % voxels activated in each of the filters
cc_mask_xpatxflt = ones(n1,n2,n3); 

if ncc == nflt+1
  cc_mask_xflt2 = ones(n1,n2,n3,npat);  % voxels activated from the concatenated time course. 
  cc_mask_xpatxflt2 = ones(n1,n2,n3);
  cc_xflt = BrikLoad(cc_fname{nflt+1});
end

% clusterize the masks
for i=1:npat
    for j=1:nflt
        % use the first line below for cc files with t values included.
      cc_mask(:,:,:,i,j) = clusterize(cc(:,:,:,(i-1)*2+1,j)>=thr|cc(:,:,:,(i-1)*2+1,j)<-thr,20,info);
     %cc_mask(:,:,:,i,j) = clusterize(cc(:,:,:,i,j)>=thr | cc(:,:,:,i,j)<-thr,20,info);
      cc_mask_xflt(:,:,:,i) = cc_mask_xflt(:,:,:,i) & cc_mask(:,:,:,i,j);
    end
    if ncc == nflt+1
    cc_mask_xflt2(:,:,:,i) = clusterize(cc_xflt(:,:,:,(i-1)*2+1)>=thr|cc_xflt(:,:,:,(i-1)*2+1)<-thr,20,info);
    end
end

for i=1:nflt
   for j=1:npat
   
     if (subbrik_mask(j))
       cc_mask_xpatxflt = cc_mask_xpatxflt & cc_mask(:,:,:,j,i);
       cc_mask_xpat(:,:,:,i) = cc_mask(:,:,:,j,i) & cc_mask_xpat(:,:,:,i);
     end

   end
end

for j=1:npat
   
     if (subbrik_mask(j) && ncc==nflt+1)
       cc_mask_xpatxflt2 = cc_mask_xpatxflt2 & cc_mask_xflt2(:,:,:,j);
     end

end

vmask = zeros(n1,n2,n3,nrois);

ts_xpat = zeros(size(data,4)/rep/npat,rep,nrois,npat,nflt); % activated across patterns
%er_xpat = zeros(size(data,4)/rep/npat,rep,nrois,npat,nflt);
ts_xflt = zeros(size(data,4)/rep/npat,rep,nrois,npat,nflt); % activated across filters
%er_xflt = zeros(size(data,4)/rep/npat,nrois,npat,nflt);
ts_all = zeros(size(data,4)/rep/npat,rep,nrois,npat,nflt);   %select all the activated voxels within each area.
%er_all = zeros(size(data,4)/rep/npat,nrois,npat,nflt);
ts_xpatxflt = zeros(size(data,4)/rep/npat,rep,nrois,npat,nflt); % activated across patterns and filters

if ncc==nflt+1
  ts_xpatxflt2 = zeros(size(data,4)/rep/npat,rep,nrois,npat,nflt); % activated across patterns and filters
  ts_xflt2 = zeros(size(data,4)/rep/npat,rep,nrois,npat,nflt); % activated across filters, as determined from concatenated time series.
end

ts_len = size(ts_xpat,1);

for i=1:nrois
    
if masktype == 0  % only include voxels map to a single roi
    umapl = (maskl(:,:,:,1)== 5);  % only map to a single roi.
    umapr = (maskr(:,:,:,1)== 5);
    if i<nrois
    vmask(:,:,:,i) = (umapl & maskl(:,:,:,2) ==i) ...
                   | (umapr & maskr(:,:,:,2)==i);
    else
     vmask(:,:,:,i) = maskl(:,:,:,1)== 2 |maskr(:,:,:,1)== 2;    % outside any roi
    end
    disp('Always use both hemispheres');
else  % include voxels mapped to multiple rois. 
    % Assign such voxels to the roi on which it has the largest number of
    % nodes.
    if strcmp(lr,'lh')
        vmask(:,:,:,i) = (maskl(:,:,:,3+nrois*2+i) == 1 );
    elseif strcmp(lr,'rh')
        vmask(:,:,:,i) = (maskr(:,:,:,3+nrois*2+i) == 1 );
    elseif strcmp(lr,'both')
     vmask(:,:,:,i) = (maskl(:,:,:,3+nrois*2+i) == 1 | maskr(:,:,:,3+nrois*2+i) == 1);
    else
        error('Unknown value for lr');
    end
end

  for ift = 1:nflt
     for ipt = 1:npat
        disp([ift ipt]);
       t_ind1 = (ipt-1)*ts_len*rep+1;
       t_ind2 = ipt*ts_len*rep;
       
       if (thr>0) % we only need to calculate this when thr > 0;
         vi = find(vmask(:,:,:,i) & cc_mask_xpat(:,:,:,ift) );
         ts_xpat(:,:,i,ipt,ift) = add_tseries(data(:,:,:,t_ind1:t_ind2,ift),vi,rep);
       
         vi = find(vmask(:,:,:,i) & cc_mask_xflt(:,:,:,ipt) );
         ts_xflt(:,:,i,ipt,ift) = add_tseries(data(:,:,:,t_ind1:t_ind2,ift),vi,rep);
       
         vi = find(vmask(:,:,:,i) & cc_mask_xpatxflt(:,:,:) );
         ts_xpatxflt(:,:,i,ipt,ift)  = add_tseries(data(:,:,:,t_ind1:t_ind2,ift),vi,rep);
        
         if ncc == nflt+1
           vi = find(vmask(:,:,:,i) & cc_mask_xflt2(:,:,:,ipt) );
           ts_xflt2(:,:,i,ipt,ift) = add_tseries(data(:,:,:,t_ind1:t_ind2,ift),vi,rep);
       
           vi = find(vmask(:,:,:,i) & cc_mask_xpatxflt2(:,:,:) );
           ts_xpatxflt2(:,:,i,ipt,ift)  = add_tseries(data(:,:,:,t_ind1:t_ind2,ift),vi,rep);
         end
       end
       
       vi = find(vmask(:,:,:,i) & cc_mask(:,:,:,ipt,ift) );
       ts_all(:,:,i,ipt,ift)  = add_tseries(data(:,:,:,t_ind1:t_ind2,ift),vi,rep);
      
       
       
     end
  end
  
end

clear data maskl maskr;

save(out_fname);

function ts = add_tseries(data,vi,rep)     

%rep:  reptitions     
     sz = [size(data,1),size(data,2),size(data,3)];
    
     duration = size(data,4)/rep;
     ts = zeros(duration,rep);     
     
     for i=1:length(vi)
    
      [i1,i2,i3]=ind2sub(sz,vi(i));
      for j=1:rep
            fst = (j-1)*duration +1;
            lst = j*duration;
            ts(:,j) = ts(:,j)+ squeeze(data(i1,i2,i3,fst:lst));
      end
     end
    
    ts = ts/length(vi);
    %er = std(ts_tmp,0,2)/sqrt(rep)/length(vi);
      
    
    
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



      
         
         
         
         


