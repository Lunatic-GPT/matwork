function sub_recon_ESPIRiT3D(prefix,step)

    lBaseResolution = readsPar(fullfile(prefix,[prefix,'.pro']),'lBaseResolution');
    lPhaseEncodingLines =readsPar(fullfile(prefix,[prefix,'.pro']),'lPhaseEncodingLines');
    lPartitions = readsPar(fullfile(prefix,[prefix,'.pro']),'lPartitions');
    lImagesPerSlab = readsPar(fullfile(prefix,[prefix,'.pro']),'lImagesPerSlab');
    lMatrix=[ lBaseResolution, lPhaseEncodingLines, lImagesPerSlab,lPartitions ];
    
   
    
    fid=fopen('ESPIRiT_script.sh','w');
    nmaps=2;
  
  
    ind=0:step:lBaseResolution;

    for i=1:length(ind)-1
        cmd=sprintf('recon_ESPIRiT3D(''%s.mat'',%d:%d,%d,[%d,%d,%d,%d])',prefix,ind(i)+1,ind(i+1),nmaps,lMatrix);
        
        fprintf(fid,'bsub -M 26 -o logESPIRiT%i.%%J matbgk "%s" logESPIRiT_mat%d\n',i,cmd,i);
    end
    
    str=num2str(ind+1);
    str=strrep(str,'  ',',');
    str=strrep(str,',,',',');
    
    fprintf(fid,'bsub -M 16 -o logESPIRiTCombine.%%J matbgk "do_combine_afterESPIRiT3D(''%s'',%d,''%s'')" logESPIRiTCombine_mat\n',prefix,nmaps,str);
    fclose(fid);
    