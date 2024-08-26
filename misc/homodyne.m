function homodyne(fmag,fph,nvox)
%homodyne(fmag,fph,nvox)
[m,info]=read_afni_sdt_images(fmag);

p=read_afni_sdt_images(fph);

z=m.*exp(1i*p);


scale=0;
mexp=zeros(4*nvox+1,4*nvox+1);


        for  i2=-2*nvox:nvox*2
            for j2=-nvox*2:nvox*2
                scale=scale+exp(-(i2^2+j2^2)/nvox/nvox);
                mexp(i2+2*nvox+1,j2+2*nvox+1)=exp(-(i2^2+j2^2)/nvox/nvox);
            end
        end
        
z2=zeros(size(z));

for k=1:size(m,3)
for i=1:size(m,1)
    for j=1:size(m,2)
      tmp=0;
        for  i2=-2*nvox:nvox*2
            for j2=-nvox*2:nvox*2
                if i+i2>size(z,1) || i+i2<1 || j+j2>size(z,2) || j+j2<1
                    continue;
                else
                    tmp=tmp+z(i+i2,j+j2,k)*mexp(i2+2*nvox+1,j2+2*nvox+1);
                end
            end
        end
        
        z2(i,j,k)=tmp;%/scale;
         
    end
    disp(i);
end
disp('');
end

prefix=strtok(fmag,'+');

WriteBrikEZ(abs(z2),info,'homodyne',[prefix,'_homo_nvox_',num2str(nvox)],'');

prefix=strtok(fph,'+');

WriteBrikEZ(angle(z2),info,'homodyne',[prefix,'_homo_nvox_',num2str(nvox)],'');

WriteBrikEZ(angle(z)-angle(z2),info,'homodyne',[prefix,'_homo_dph_nvox_',num2str(nvox)],'');

