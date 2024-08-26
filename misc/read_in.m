function img_dat = read_in(filename,NP,NF,tot_slice,slice_idx,time_pt,num_offset)
%read in afni dataset explicitly using fread() instead of using BrikLoad() 
        header_size=0; %afni_orig.BRIK
        dat_size=2*NF*NP;
        img_dat=zeros(NP,NF,time_pt);
	for k=1:time_pt
         offset=header_size+dat_size*num_offset+dat_size*(slice_idx-1)+dat_size*(k-1)*tot_slice;
         fid=fopen(filename,'r');
         status=fseek(fid,offset,'bof');
         [raw_dat,COUNT]=fread(fid,dat_size,'short');
	 fclose(fid);
         if COUNT ~= dat_size
          sprintf('error in reading');
         end
         %image construction
         for j=1:NP %phase 
          for i=1:NF %frequency
           img_dat(j,i,k)=raw_dat(i+(j-1)*NP);
          end
         end
	end