function trace2dset(fname)

 f=fopen(fname,'r');
 
 a=textscan(f,'%s');
 fclose(f);
 
 
 for i=1:length(a{1})
     if strcmp(a{1}{i},'<samplespacing')  
         for j=1:3
         tmp=a{1}{i+j};
         [tmp,rem]=strtok(tmp,'=');
         ss(j)=str2num(rem(3:end-1));
         end
     elseif strcmp(a{1}{i},'<imagesize')
          for j=1:3
           tmp=a{1}{i+j};
           [tmp,rem]=strtok(tmp,'=');
           if j<3
            sz(j)=str2num(rem(3:end-1));
           else
             sz(j)=str2num(rem(3:end-3));
           end
          end
          data=zeros(sz);
     elseif strcmp(a{1}{i},'<point')
         for j=1:3
           tmp=a{1}{i+j};
           [tmp,rem]=strtok(tmp,'=');
            ind(j)=str2num(rem(3:end-1));
         end
         data(ind(1),ind(2),ind(3))=1;
     end
     
     
 end
 
 ss([2,1])=ss([1,2]);
 data=permute(data,[2,1,3]);
 data=flip(data,1);
 %data=flip(data,2);
 
 info.DELTA=ss;
  %  info.DELTA. voxel size (default [1,1,1]);  negative if
  %  different from [0 3 4]
  %  info.ORIENT_SPECIFIC(default [0 3 4]);(R2L,A2P,I2S; afni convention);  
  fname=strtok(fname,'.');
 write_afni(data,fname,info);

 