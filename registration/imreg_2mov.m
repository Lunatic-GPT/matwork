function [res,res2]=imreg_2mov(fix,mov,mov2,do_circshift,out_name)
% [res,res2]=imreg_2mov(fix,mov,mov2,do_circshift,out_name)
% use the transformation of mov ->fix to transform mov2 as well
% res: registered image for mov; 
% res2: reigstered image for mov2
fix=ri(fix);
mov=ri(mov);
mov2=ri(mov2);

if do_circshift
    [d,output]=imreg_dft(fix(:,:,:,1),mov(:,:,:,1));
    res=circshift(mov2,[round(output(3:4)),0,0]);
    disp(output);
else
 [res,xfm]=imreg_matlab(fix(:,:,:,1),mov(:,:,:,1));
if ~isempty(mov2) 
 res2=imwarp(uint16(mov2),xfm,'nearest','OutputView',imref2d(size(mov2)));
end

end

if exist('out_name','var')
save(out_name,'res','res2');
end
