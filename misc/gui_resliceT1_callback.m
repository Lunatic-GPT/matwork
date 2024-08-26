function gui_resliceT1_callback(p)

T1=get_fpattern(p,'T1 image');
T1=str2cell(T1);
for i=1:length(T1)

    if ~isempty(strfind(T1{i},'.'))
        continue;
    end
    d=ri(T1{i},[]);

o_voxsize=dcmDimCenter(T1{i});

u_voxsize=get(p,'new dimension (mm)');
udsz=get(p,'new matrix');
center_u2o=get(p,'center shift');

if(get(p,'flip xy'))
    d=permute(d,[2,1,3]);
end


if(get(p,'flip ap'))
    d=flip(d,2);
end


if(get(p,'flip lr'))
    d=flip(d,1);
end

for j=1:size(d,3)
 [v2(:,:,j),sl_u2o]=reslice_extend(d(:,:,j),o_voxsize,u_voxsize,udsz,center_u2o,true);
end

fname=filename(T1{i});

save([fname,'_reslice'],'v2');
disp('reslice done');
end