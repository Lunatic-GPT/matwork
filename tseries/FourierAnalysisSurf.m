function brikdata = FourierAnalysisSurf(in_fname,num_img_pd,out_prefix)
%
% FourierAnalysis(in_fname,clip_fname,num_img_pd[,out_prefix])
%
% 
% in_fname:     1D dataset name. Use column 0 as node index column.
% clip_fname: 1D file for clipvalue and average
% num_img_pd: number of TR's in each period.
% out_prefix:   prefix of the 1D functional data;
%               do not save if leave blank.
% brikdata: 1. phase 2. corr coef 3. perc modu
%
% 
% 2-9-2009, XZ	
	tic;

	%start data analysis
  
  	%Read in func data from r2+orig.BRIK 
       % filename='r60dreg_t+orig.BRIK';
       % filename='r60d23inreg_t+orig.BRIK';
        %filename='w30dreg+orig.BRIK';
         %filename = 'w30dtshift+orig.BRIK';
        func_dat=textread(in_fname,'','commentstyle','shell');
        nn = size(func_dat,1); %number of nodes.
      
       sti_ir=cal_sti_ir(func_dat(:,2:end),num_img_pd);
 
        
	%calculate the stimulus phase.  The angle increases in the counter-clockwise
	%direction
	   
       cc_mag = sqrt(sti_ir(:,1).^2+sti_ir(:,2).^2); 
       
       
    %save data out
       
       if exist('out_prefix','var')
      
         fid = fopen([out_prefix,'.1D.dset'],'w');
         fprintf(fid,'#index\tphase\tcorr coef\n');
         for i=1:nn
         fprintf(fid,'%d\t%f\t%f\n',func_dat(i,1),sti_ir(i,3),cc_mag(i));
         end

         fclose(fid);
           
       end
       disp(sprintf('Program total running time: %f s',toc)); 
	
	
 
     
% calculate the Fourier transformation
function sti_ir=cal_sti_ir(func_dat,num_img_pd)
	
    nt = size(func_dat,2);
        
    if (mod(nt,num_img_pd) ~= 0)
	 disp('error nt or num_img_pd');
	else
	 num_cyc=floor(nt/num_img_pd);
    end
    
    nn=size(func_dat,1);
	sti_ir=zeros(nn,4);

    for i=1:nn
	       
             tmp=squeeze(func_dat(i,:));
              
             if isempty(find(tmp~=0,1))                        
	            continue;
             end
             norm = (nt-1)*sqrt(nt/(nt-1)/2)*std(tmp);
             if norm == 0
                 continue;
             end
	         ft=fft(tmp);
	         sti_ir(i,1)=imag(ft(num_cyc+1))/norm;
	         sti_ir(i,2)=real(ft(num_cyc+1))/norm;    
             
              tmpAngle = -angle(ft(num_cyc+1))*180/pi; 
	          sti_ir(i,3)= mod(tmpAngle,360);
              sti_ir(i,4)= abs(ft(num_cyc+1));

    end
    
    