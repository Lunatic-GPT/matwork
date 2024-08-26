function extract_diffusionvec_from_dcm(dcm,outname)

extp(dcm);
pro=[strtok2(dcm,'.'),'.pro'];
counter=1;
vec_len=readsPar(pro,'sDiffusion.sFreeDiffusionData.lDiffDirections');
v=zeros(vec_len-1,3);

for i=1:vec_len-1
    sag=sprintf('sDiffusion.sFreeDiffusionData.asDiffDirVector[%d].dSag',i);
    cor=sprintf('sDiffusion.sFreeDiffusionData.asDiffDirVector[%d].dCor',i);
    tra=sprintf('sDiffusion.sFreeDiffusionData.asDiffDirVector[%d].dTra',i);
     x=readsPar(pro,sag);
   y=readsPar(pro,cor);
   z=readsPar(pro,tra);

    if exist('vec_len','var') && counter>vec_len
        break;
    end
      
   v(i,:)=[x,y,z];
   
   disp([i,v(i,:)]);

   
       
end
disp(sos(v,2));
write_DiffusionVectors(v,outname);