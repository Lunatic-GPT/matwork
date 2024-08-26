function show_cd(data,roi,range,range_angle)

if ~exist('roi','var')
    roi=ones(size(data));
end

if ~exist('range','var')
  range=min_max([imag(data(:));real(data(:))]);
end

if ~exist('range_angle','var')
    range_angle=[-180,180];
end

range_angle=range_angle/180*pi;

data(roi==0)=0;


 an_data=(angle(data)-range_angle(1))/(range_angle(2)-range_angle(1))*(range(2)-range(1))+range(1);
 

  im=cat(3,real(data),imag(data),an_data);
  
  im=repmat2(im,20);
  if ~exist('range','var')
      range=min_max(im(:));
  end
  imshow4(im,range,[1,3],2);
  