function [reg,tform]=imreg_matlab(fix,mov,fout,threeD,do_reg,type)
%reg=imreg_matlab(fix,mov,fout[,threeD])
% threeD: default false
% do_reg: if false, only tform will be calculated; reg will be empty.
% type: 'translation','rigid', 'similarity','affine'; default: rigid
% tform is the transform matrix that can be used to align other images using the same tform.
% Example: imwarp(mov,tform,'OutputView',imref2d(size(fix)));


if isa(fix,'char')
    fix=ri(fix);
end
fix=fix(:,:,:,1);
if isa(mov,'char')
    mov=ri(mov);
end

if ~exist('threeD','var')
    threeD=false;
end

if ~exist('do_reg','var')
    do_reg=true;
end

if ~exist('type','var')
    type='rigid';
end

[optimizer, metric] = imregconfig('multimodal');

% Tune the properties of the optimizer to get the problem to converge
% on a global maxima and to allow for more iterations.
%optimizer.InitialRadius = 0.0005;  % this parameter is important for registeration accuracy.
%optimizer.Epsilon = 1.5e-5;
%optimizer.GrowthFactor = 1.01;
%optimizer.MaximumIterations = 600;

% Align the moving image with the fixed image

out=zeros(size(mov));
tic;
for i=1:size(mov,4)

    if threeD
        trans_mat_s = imregtform(mov(:,:,:,i),fix, type, optimizer, metric);


        if do_reg
            out(:,:,:,i)=imwarp(mov(:,:,:,i),tform(i),'OutputView',imref2d(size(fix)));
        end
    else
        for j=1:size(mov,3)
            trans_mat_s = imregtform(mov(:,:,j,i),fix, type, optimizer, metric);
            if do_reg
                out(:,:,j,i)=imwarp(mov(:,:,j,i),tform(i),'OutputView',imref2d(size(fix)));
            end
        end
    end

    trans_mat = trans_mat_s.T;

    trans_mat = trans_mat';
    tran(i,:)=trans_mat(1:3,4)';
    mat = [trans_mat(1) trans_mat(2) trans_mat(3);
        trans_mat(5) trans_mat(6) trans_mat(7);
        trans_mat(9) trans_mat(10) trans_mat(11)];
    angle(i,:) = get_rot_axis(mat);



    time_left(i,size(mov,4),toc);
end

if do_reg
    reg=abs(out);
else
    reg=[];
end


if exist('fout','var')&&~isempty(fout)
    save(fout,'reg','tran','angle');
end

%write_afni(abs(out),fout);