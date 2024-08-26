function image_contour(varargin)

%fil = ones(3,3);
img = varargin{1};
img = floor(img*100/max(img(:)));
img(img<=2) = 2;

cm=gray(100);
cm(1,:)=[1,0,0];
colormap(cm);
for i=2:nargin
     %roi_morph = bwmorph(varargin{i},'remove');
     %TG1 = filter2(fil,roi_morph);
     %TG1(TG1>0)=1;  
    
     %gap_BW = imfill(TG1,'holes')-varargin{i};
      gap_BW=bwmorph(varargin{i},'remove');
     img(gap_BW>0) = 1;     
     
end

image(img);
axis image off;