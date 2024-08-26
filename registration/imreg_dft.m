function [reg,output]=imreg_dft(fix,mov,fout)
%imreg_dft(fix,mov,fout)

if isa(fix,'char')
fix=BrikLoadf(fix);
end
if isa(mov,'char')
mov=BrikLoadf(mov);
end

out=zeros(size(mov));

for i=1:size(mov,4)
    for j=1:size(mov,3)
      [output, reg] = dftregistration(fft2(fix(:,:,j)),fft2(mov(:,:,j,i)),20);
      out(:,:,j,i)=ifft2(reg);
    end 
end

reg=abs(out);
if exist('fout','var')
save(fout,'reg');
end
%write_afni(abs(out),fout);