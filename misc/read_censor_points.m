function plist=read_censor_points(fname)
% the indices in fname are 0 based, to be consistent with afni convention.
fid = fopen(fname,'r');
s = textscan(fid,'%s');
plist = s{1};
for i=1:length(plist)
    plist{i} = str2num(plist{i});
    if any(plist{i}<0)
        plist{i} = [];
    end
end