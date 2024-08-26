function brikdata = FourierAnalysis(in_fname,clip_fname,clipLevel, num_img_pd,out_prefix)
%
% FourierAnalysis(in_fname,clip_fname,clipLevel,num_img_pd[,out_prefix])
%
% 
% in_fname:     afni dataset name
% clip_value: data set for calculate the average signal intensity
% num_img_pd: number of TR's in each period.
% out_prefix:   prefix of the output functional data, (phase and corrcoef);
%               do not save if leave blank.
% brikdata: 1. phase 2. corr coef 3. perc modu
%
% 1-8-2009, XZ	
	tic;

	%start data analysis
  
  	%Read in func data from r2+orig.BRIK 
       % filename='r60dreg_t+orig.BRIK';
       % filename='r60d23inreg_t+orig.BRIK';
        %filename='w30dreg+orig.BRIK';
         %filename = 'w30dtshift+orig.BRIK';
        [func_dat,info] = BrikLoad(in_fname);
        nf = size(func_dat,1);
        np = size(func_dat,2);
        num_loc = size(func_dat,3);
        nt = size(func_dat,4);
             
      
       sti_ir=cal_sti_ir(func_dat,num_img_pd);
 
        
	%calculate the stimulus phase.  The angle increases in the counter-clockwise
	%direction
	   
       cc_mag = sqrt(sti_ir(:,:,:,1).^2+sti_ir(:,:,:,2).^2); 
        
       [clip_dat,info] = BrikLoad(clip_fname);
        ave_dat = mean(clip_dat,4);
        
        pm = zeros(nf,np,num_loc);
        % calculate the peak to peak percent modulation 
        for i=1:nf
            for j=1:np
                for k=1:num_loc
                    if ave_dat(i,j,k)>=clipLevel
                        pm(i,j,k) = sti_ir(i,j,k,4)*4/nt/ave_dat(i,j,k);
                    end
                end
            end
        end
        
        
    %save data out
    
       brikdata = zeros(nf,np,num_loc,3);
       brikdata(:,:,:,1) = sti_ir(:,:,:,3);
       brikdata(:,:,:,2) = cc_mag;
       brikdata(:,:,:,3) = pm;
       
       if exist('out_prefix','var')
       info.DATASET_RANK(2) = 3;  % the number of subbricks to save.
       info.BRICK_TYPES=[3 3 3];  %1: short, 3: float
       info.BRICK_LABS = 'phase~ corr coef~perc modu~';
       info.IDCODE_STRING = 'ring wedge maps';
       
       info.TYPESTRING = '3DIM_HEAD_FUNC';
       info.SCENE_DATA = [0 0 1];  % first 0: orig view; second 0: 1 value; 1: matches TYPESTRING 
       
       info.BRICK_STATS=[];  %minimum and maximum values of the subbrick;
       
       info.BRICK_FLOAT_FACS = [];
       
       info.IDCODE_DATE = date;
       info.TAXIS_NUMS=[];
       
       %opt.Prefix = 'retino_map_phase2+orig';
       opt.Prefix = out_prefix;
       opt.OverWrite = 'y';
       WriteBrik(brikdata,info,opt);
       end
       disp(sprintf('Program total running time: %f s',toc)); 
	
	
    % generate reference function and return the reference functions
	function ref_dat = read_ref(fnames)
        for i=1:length(fnames)
            tmp = textread(fnames{i});
            ref_dat(:,i)= tmp/max(tmp);
        end
	
	
     
% calculate the Fourier transformation
function sti_ir=cal_sti_ir(func_dat,num_img_pd)
	
    nt = size(func_dat,4);
        
    if (mod(nt,num_img_pd) ~= 0)
	 disp('error nt or num_img_pd');
	else
	 num_cyc=floor(nt/num_img_pd);
    end
    
    nf=size(func_dat,1);
    np=size(func_dat,2);
    num_loc=size(func_dat,3);
   
	sti_ir=zeros(nf,np,num_loc,4);

    for i=1:nf
	   for j=1:np
         for k=1:num_loc
	       
             tmp=squeeze(func_dat(i,j,k,:));
              
             if isempty(find(tmp~=0,1))                        
	            continue;
             end
             norm = (nt-1)*sqrt(nt/(nt-1)/2)*std(tmp);
             if norm == 0
                 continue;
             end
	         ft=fft(tmp);
	         sti_ir(i,j,k,1)=imag(ft(num_cyc+1))/norm;
	         sti_ir(i,j,k,2)=real(ft(num_cyc+1))/norm;    
             
              tmpAngle = -angle(ft(num_cyc+1))*180/pi; 
	          sti_ir(i,j,k,3)= mod(tmpAngle,360);
              sti_ir(i,j,k,4)= abs(ft(num_cyc+1));
         end
       end
    end
    
    