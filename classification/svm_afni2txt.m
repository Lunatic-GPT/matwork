function svm_afni2txt(fname,fmask,category,txtname,frames)
% function svm_afni2txt(fname,fmask,category,txtname,frames)
% fname: name of the afni dataset
% fmask: mask in afni format.  Voxels with non-zero values in the mask will be used in the patterns.
% category: category labels of the briks
% txtname: file name of the output text file.
% frames: used for subbrik and category label selection. 1 based.  use []
% to select all frames.
%NOTE: the category label length and the total number of subbriks in fname
%should be the same.  


labels = load(category);
[err,info] = BrikInfo(fname);
if info.DATASET_RANK(2) ~=length(labels)
    error('Category length error');
end

if exist('frames','var') && ~isempty(frames);
    Opt.Frames = frames;
    d =BrikLoad(fname,Opt);
    labels = labels(frames);
else
    d = BrikLoad(fname);
end

mask = BrikLoad(fmask);
ind = find(mask>0);
% generate training file.
sz = size(d);

fid = fopen(txtname,'w');

  t_tr = find(labels<9999);
  train_str = char(zeros(1,100+19*length(ind)*length(t_tr)));
  p = 0;
  fprintf('Reading %d examples ', length(t_tr));
for t=1:length(t_tr)
    
  str = sprintf('%d',labels(t_tr(t)));  
  train_str(p+1:p+length(str)) = str;
  p = p +length(str);
  for j=1:length(ind)
      [i1,i2,i3] = ind2sub(sz(1:3),ind(j));
      str = sprintf(' %d:%f',j,d(i1,i2,i3,t_tr(t)));
      train_str(p+1:p+length(str)) = str;
      p = p +length(str);
  end
  train_str(p+1) = sprintf('\n');
  p = p+1;
  fprintf('.');
end
fprintf('\n');
fprintf(fid,'%s',train_str);
fclose(fid);
