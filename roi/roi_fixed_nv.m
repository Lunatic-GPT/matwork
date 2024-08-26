function [m,thr] = roi_fixed_nv(varargin)
% [m,thr] = roi_fixed_nv(cc1[,cc2,cc3...],brainmask,n)
% or [m,thr] = roi_fixed_nv(cc1[,cc2,cc3...],brainmask,n,prefix)
% the program will determine the threshold for cc value such that only n
% voxels are included in the final mask m.  m cosists of voxels passing
% thresholding in all cc data.
% 
% cc1 - ccn:  mask data or an afni brik 
% n: total number of voxels to include in the final mask. 
% if use the second form, cc might be file names.

if isa(varargin{end},'char')
    nv = varargin{end-1};
    ncc = nargin-3;
    brainmask = varargin{end-2};
else
    nv = varargin{end};
    ncc = nargin-2;
    brainmask = varargin{end-1};
end

bm = BrikLoad(brainmask);
bm = (bm>0);

for i=1:ncc
    if isa(varargin{i},'char')
        [cc{i},info] = BrikLoadf(varargin{i});
        cc{i}(~bm) = 0;
    else
        cc{i} = varargin{i};
        cc{i}(~bm) = 0;
    end
end


for i=1:ncc
  [cc_sort{i},ind{i}] = sort(cc{i}(:),'descend');
end



n = nv;
while 1
cmmn = ind{1}(1:n);
for i=2:ncc
    cmmn = intersect(cmmn,ind{i}(1:n));
end

if length(cmmn(:))>=nv
    break;
end

n=n+1;
end


sz = size(cc{1});
m = false(sz(1:3));

for i=1:nv
    [i1,i2,i3] = ind2sub(sz(1:3),cmmn(i));
    m(i1,i2,i3) = true;
end

thr = cc_sort{1}(n);  % initial threshold
for i=2:ncc
    if thr>cc_sort{i}(n)
        thr = cc_sort{i}(n);
    end
end

    

if isa(varargin{end},'char')
history =   'roi_fixed_nv(';
for i=2:nargin-2
    history = [history,varargin{i},','];
end

history = sprintf('%s%d,%s)',history,varargin{end-1},varargin{end});


WriteBrikEZ(m,info,history,varargin{end});
end

   



    
    
