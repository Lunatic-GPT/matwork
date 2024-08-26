function [img, dimStr, numDimStr] = rdSdt( fileNameRoot ,endian)
%[img, dimStr, numDimStr] = rdSdt( fileNameRoot[,endian] )
format compact

sprPath = [fileNameRoot, '.spr'];
fp = fopen( sprPath ); 
done=0;
while( done == 0 )
   line = fgetl(fp);
   if (line == -1)
      done = 1;
   else
      [attr, val] = strtok(line,':');
      switch attr
         case 'numDim'
            numDimStr = strtok(val,':');
         case 'dim'
            dimStr = strtok(val,':');
         case 'dataType'
            dataTypeStr = strtok(val,':');
         case 'Real2WordScale'
            Real2WordScale = strtok(val,':');
          case 'endian'
              endian = strtrim(strtok(val,':'));
      
      end
   end
end
fclose( fp );


if ~exist('endian','var')
    endian = 'ieee-le';
end

switch deblank(strjust( dataTypeStr ,'left'))
   case 'BYTE'
		precision = 'uchar';
   case 'REAL'
		precision = 'float32';
        sdtPath = [fileNameRoot, '.sdt'];
        fp = fopen( sdtPath ,'r',endian); 
        dimArr = str2num(dimStr);
        
        img = fread( fp, prod(dimArr), precision);
        if size(dimArr,2) == 4
            img = reshape( img, dimArr(1),dimArr(2),dimArr(3),dimArr(4) );
        else if size(dimArr,2) == 3
            img = reshape( img, dimArr(1),dimArr(2),dimArr(3) );
            else
            img = reshape( img, dimArr(1),dimArr(2) );
            end
        end
  
        
   case 'WORD'
		precision = 'int16';
      %  scale = str2num(Real2WordScale);
      scale =1; 
      sdtPath = [fileNameRoot, '.sdt']
        fp = fopen( sdtPath ); 
        dimArr = str2num(dimStr);
        img = zeros( dimArr );
        img = fread( fp, prod(dimArr), precision,endian);
        if size(dimArr,2) == 4
            img = reshape( img, dimArr(1),dimArr(2),dimArr(3),dimArr(4) );
        else if size(dimArr,2) == 3
            img = reshape( img, dimArr(1),dimArr(2),dimArr(3) );
            else
            img = reshape( img, dimArr(1),dimArr(2) );
            end
        end
        %img = reshape( img, dimArr(1),dimArr(2),dimArr(3),dimArr(4) );
        img = img/scale;
   otherwise
		error( 'Can not handle dataType')
end
%img=flipdim(img,2);
fclose( fp );

