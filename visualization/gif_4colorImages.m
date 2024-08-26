function gif_4colorImages(d,map,prefix,delayTime,crop)
%  gif_4colorImages(d,map,prefix,delayTime)
% data in gif file are saved as uint8, although the input data can be in other formats

%d: can be 3D matrix or a cell of tiff file names.  In the latter case, map is
%empty.
%map values should be between 0 and 1.

if ~exist('delayTime','var')
    delayTime = 0.1;
end


if isa(d,'cell')
    fname=d;
    for i=1:length(fname)
        if i==1
            d=imread(fname{i});
        else        
            d(:,:,:,i)=imread(fname{i});
        end
    end
    [d,map]=convert2Map(d);
    map=double(map)/255;
end

d=squeeze(d);
if exist('crop','var')
    lt=crop(1);
    rt=crop(2);
    up=crop(3);
    dn=crop(4);
end

if exist('crop','var')
    d=d(up+1:end-dn,lt+1:end-rt,:);
end

for n=1:size(d,3)
    if n == 1
        imwrite(double(d(:,:,n)),map,[prefix,'.gif'],'gif','LoopCount',Inf,'DelayTime',delayTime);
    else
        imwrite(double(d(:,:,n)),map,[prefix,'.gif'],'gif','WriteMode','append','DelayTime',delayTime);
    end
end

function [ind,map]=convert2Map(d)

% d is row*col*3*frame

d2=permute(d,[1,2,4,3]);

sz=size(d2);

d2=reshape(d2,[prod(sz(1:3)),3]);

[map,~,ind]=unique(d2,'rows');

ind=reshape(ind,sz(1:3));














