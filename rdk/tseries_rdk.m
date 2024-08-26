function tseries_strpsq(p,thr,lr,masktype)
% 
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
% 5-26-2009: added input argument ntrl_ts, so that only the last ntrl_ts trials
%            are used for time series.  The remaining trials were
%            supposedly used in correlation coefficient calculation.
% 9-16-2009: time series of each pattern is stored in a seperate file.  Before, all
% patterns of the same filter were stored in the same file.  Removed npat
% as an input variable.
% old format:
% tseries(ts_fname,cc_fname,out_fname,thr,rep,npat,[ntrl_ts,lr,subbrik_mask,masktype,vml_name,vmr_name])


tic;

out_fname = get(p,'output file');


%
pm_strpsq = get(p,'pattern mask for cc threshold (strpsqr)');
pm_strp = get(p,'pattern mask for cc threshold (stripes)');
pm_strpsq=str2num(pm_strpsq);
pm_strp = str2num(pm_strp);

ts_fname = get(p,'ordered time series');
ts_fname = str2list(ts_fname);
npat = length(ts_fname);
rep = get(p,'repititions per stimulus');

cc_fname = get(p,'cc');


vml_name = get(p,'left rois');
vmr_name = get(p,'right rois');
%lr = get(p,'which roi to use');

if ~exist('masktype','var')
    masktype = 1;
end

if ~exist('subbrik_mask','var')
    subbrik_mask = ones(1,npat);
end


[vmask,rois_l] = get_vmask(vml_name,vmr_name,lr,masktype);
nrois = size(vmask,4);

n1 = size(vmask,1);
n2 = size(vmask,2);
n3 = size(vmask,3);

    for ipat = 1:npat
      if ipat==1
        data = BrikLoad(ts_fname{ipat});
      else
          data_temp2 = BrikLoad(ts_fname{ipat});
          data = cat(4,data,data_temp2);
      end  
    end
    
cc_mask = zeros(n1,n2,n3,npat);
cc_mask_xpat = ones(n1,n2,n3);
[cc,info] = BrikLoad(cc_fname);
% clusterize the masks
for i=1:npat
        % use the first line below for cc files with t values included.
      cc_mask(:,:,:,i) = clusterize(cc(:,:,:,(i-1)*2+1)>=thr|cc(:,:,:,(i-1)*2+1)<-thr,20,info);
end

   for j=1:npat
     if (pm_strpsq(j))
       cc_mask_xpat(:,:,:) = cc_mask(:,:,:,j) & cc_mask_xpat(:,:,:);
     end
   end

cc_fname2{1} = get(p,'cc (gray)');

ts_xpat = zeros(size(data,4)/rep/npat,rep,nrois,npat); % activated across patterns
ts_all = zeros(size(data,4)/rep/npat,rep,nrois,npat);   %select all the activated voxels within each area.

ts_len = size(ts_xpat,1);

for i=1:nrois
    
     for ipt = 1:npat
        fprintf('rois %d/%d, pattern %d/%d\n',i,nrois,ipt,npat);
       t_ind1 = ipt*ts_len*rep-ts_len*rep+1;
       t_ind2 = ipt*ts_len*rep;
       
       if (thr>0) % we only need to calculate this when thr > 0;
         
         ts_xpat(:,:,i,ipt) = add_tseries(data(:,:,:,t_ind1:t_ind2),vmask(:,:,:,i) & cc_mask_xpat,rep);
         ts_xflt_strp(:,:,i,ipt) = add_tseries(data(:,:,:,t_ind1:t_ind2),vmask(:,:,:,i) & cm_xflt_strp(:,:,:,mod(ipt-1,3)+1),rep);
         ts_xflt2_strp(:,:,i,ipt) = add_tseries(data(:,:,:,t_ind1:t_ind2),vmask(:,:,:,i) & cm_xflt2_strp(:,:,:,mod(ipt-1,3)+1),rep);
         ts_xpatxflt_strp(:,:,i,ipt) = add_tseries(data(:,:,:,t_ind1:t_ind2),vmask(:,:,:,i) & cm_xpatxflt_strp,rep);
         ts_xpatxflt2_strp(:,:,i,ipt) = add_tseries(data(:,:,:,t_ind1:t_ind2),vmask(:,:,:,i) & cm_xpatxflt2_strp,rep);
        end
       
       ts_all(:,:,i,ipt)  = add_tseries(data(:,:,:,t_ind1:t_ind2),vmask(:,:,:,i) & cc_mask(:,:,:,ipt) ,rep);
      
     end
  
end

clear data maskl maskr data_temp data_temp2;

save(out_fname);
disp([mfilename ' finish in ', num2str(toc), ' s']);
       


function ts = add_tseries(data,vi,rep)     

%rep:  reptitions     
     
     duration = size(data,4)/rep;   
    
     data = data.*repmat(vi,[1,1,1,size(data,4)]);
     data = sum(sum(sum(data,1),2),3);
     data = squeeze(data);
     ts = reshape(data,duration,rep);
    ts = ts/length(find(vi));
    %er = std(ts_tmp,0,2)/sqrt(rep)/length(vi);
      
    



      
         
         
         
         


