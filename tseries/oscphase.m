function oscphase(fname,num_cyc)
% oscphase(fname,num_cyc)
	tic;
        
        [func_dat,info] = BrikLoad([fname,'+orig']);
        sz = size(func_dat);
        sti_ir_wdg = zeros([sz(1:3),3]);
        
       
        for i=1:size(func_dat,1)
            for j=1:size(func_dat,2)
                for k=1:size(func_dat,3)
                 sti_ir_wdg(i,j,k,:)=cal_sti_ir(func_dat(i,j,k,:),num_cyc);
                end
            end
        end
  
         
	%calculate the stimulus phase.  The angle increases in the counter-clockwise
	%direction
	   
       cc_mag_wdg = sqrt(sti_ir_wdg(:,:,:,1).^2+sti_ir_wdg(:,:,:,2).^2); 
    
       brikdata = zeros([sz(1:3),2]);
       brikdata(:,:,:,1) = sti_ir_wdg(:,:,:,3);
       brikdata(:,:,:,2) = cc_mag_wdg;
        
        WriteBrikEZ(brikdata,info,'',[fname,'_oscphase'],'phase~cc~');
       disp(sprintf('Program total running time: %f s',toc)); 
	

     
% calculate correlation coefficient with fft to speed up the code.
%calculate correlation coefficients with sin and cos, and the phase of the stimulus
%activation time course.
function sti_ir=cal_sti_ir(func_dat,num_cyc)
	
    nt = length(func_dat);
        
    
	          sti_ir=zeros(3,1);

  
             tmp=squeeze(func_dat(:));
              
             if isempty(find(tmp~=0,1))                        
	            return;
             end
             norm = (nt-1)*sqrt(nt/(nt-1)/2)*std(tmp);
             if norm == 0
                 return;
             end
	         ft=fft(tmp);
	         sti_ir(1)=imag(ft(num_cyc+1))/norm;
	         sti_ir(2)=real(ft(num_cyc+1))/norm;    
             tmpAngle = -angle(ft(num_cyc+1))*180/pi;
             sti_ir(3)=mod(tmpAngle,360);
              
 

