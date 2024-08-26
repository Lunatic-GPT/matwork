function ccmaskROIs(mask,lr,vml_name,vmr_name)
% ccmaskROIs(mask[,lr,vml_name,vmr_name])
% ts_fname:  data sets of the time series
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

tic;
if ~exist('vml_name','var')
    vml_name = 'brainmask_sroi2v_lh+orig';
end

if ~exist('vmr_name','var')
    vmr_name = 'brainmask_sroi2v_rh+orig';
end

if ~exist('lr','var')
    lr = 'both';
end

[maskl,info] = BrikLoad(vml_name);
rois_l = roi_search(info.BRICK_LABS);
[maskr, info] = BrikLoad(vmr_name);
rois_r = roi_search(info.BRICK_LABS);

nrois = length(rois_l);

if (nrois ~=length(rois_r))
    error('The left and rigth hemispheres have different rois');
end


for i=1:nrois
    sub_sel = sprintf('-a %s',mask); 
    switch lr
        case 'lh'
            sub_sel = sprintf('%s -b %s''[%d]''',sub_sel,vml_name,3+nrois*2+i-1);
            expr = 'notzero(a)*equals(b,1)';
      %  vmask(:,:,:,i) = (maskl(:,:,:,3+nrois*2+i) == 1 );
        case 'rh'
         sub_sel = sprintf('%s -b %s''[%d]''',sub_sel,vmr_name,3+nrois*2+i-1);
         expr = 'notzero(a)*equals(b,1)';
        case 'both' 
         % vmask(:,:,:,i) = (maskl(:,:,:,3+nrois*2+i) == 1 | maskr(:,:,:,3+nrois*2+i) == 1);
          sub_sel = sprintf('%s -b %s''[%d]'' -c %s''[%d]''',sub_sel,vml_name,3+nrois*2+i-1,vmr_name,3+nrois*2+i-1);
          expr = 'notzero(a)*ispositive(equals(b,1)+equals(c,1))';
        otherwise 
        error('Unknown value for lr');
    end
    prefix = strtok(mask,'+');
    cmd = sprintf('3dcalc -prefix %s_%s%s %s -expr ''%s''',prefix,rois_l{i},lr,sub_sel,expr);
    unix(cmd);
end

disp([mfilename ' finish in ', num2str(toc), ' s']);
       



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
      
         
         
         
         


