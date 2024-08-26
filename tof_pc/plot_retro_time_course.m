function [pulsa_all,vmean_all]=plot_retro_time_course(v_all,v_all_same,v_all_nobg,h1,h2,nshift,ev,ind_min,ind_max)
%plot_retro_time_course(v_all,v_all_same,v_all_nobg,h1,h2,nshift,ev)
% file: phase file or mag file in units of degree
% if mag file: venc should be 'tof', otherwise the venc value
nt=size(v_all,2);
nshift=mod(nshift,nt);
if ~isempty(h1)
    figure(h1);
    hold off;
    
    nv=size(v_all,1);
    ncol=ceil(sqrt(nv));
    
    nrow=ceil(nv/ncol);
    
    ref=mean(v_all,1);
    ref=ref/mean(ref);
    clr='rgb';
    
    
    for i=1:nv
        subplot(nrow,ncol,i);
        hold off;
        for j=1:3
            if j==1
                if isempty(v_all)
                    continue;
                end 
                vtmp=v_all(i,:);
                 tmp=v_all(i,:);
                
            elseif j==2
                if isempty(v_all_same)
                    continue;
                end
                vtmp=v_all_same(i,:);
                 tmp=[tmp,vtmp];
                
            else
                if isempty(v_all_nobg)
                    continue;
                end
               
                
                vtmp=v_all_nobg(i,:);
                 tmp=[tmp,vtmp];
            end
            
            vtmp=circshift(vtmp,[0,nshift]);
            
            plot([vtmp,vtmp(1)],clr(j));
            
            hold on;
        end
        
        ref_tmp=circshift(ref,[0,nshift])*mean(v_all(i,:));
        plot([ref_tmp,ref_tmp(1)],'k');
        
        %  legend('ind','same','no','mean');
        
       
        ylim(min_max(tmp)+[-0.1,0.1]);
        
        xlim([0.8,size(v_all,2)+1+0.2]);
        
        yl=ylim;
        
        plot((nshift+1)*[1,1]/nt,[yl(1),yl(1)+0.1],'r-');
        
    end
end

tlt={'indiv bg','same bg','no bg corr'};

        pulsa_all=zeros(1,2);
vmean_all=zeros(1,2);
same_panel=true;
figure(h2);

        fid=fopen('plot_retro_time_course_output.txt','w');
for j=1:2  % normalized vs un-normalized
    
    if same_panel
        subplot(1,2,j);
        hold off;
    else
        subplot(3,2,i*2-2+j);
    end

    hplot=[];
    lgn={};
    for i=1:2
        
        if i==1
            v_all_tmp=v_all';
            lgn{end+1}='indv bg';
        elseif i==2
            v_all_tmp=v_all_same';
            
            lgn{end+1}='same bg';
        else
            v_all_tmp=v_all_nobg';
            
            lgn{end+1}='no bg';
        end
        if isempty(v_all_tmp)
            continue;
        end
           ind=any(isnan(v_all_tmp),1);
        v_all_tmp(:,ind)=[];
        v_all_tmp=circshift(v_all_tmp,[nshift,0]);
        
        nt=size(v_all_tmp,1);
        nv=size(v_all_tmp,2);
        if j==1
             y=mean(v_all_tmp,2);
          %   ey=std(v_all_tmp,[],2)/sqrt(size(v_all_tmp,2));  %calculate
          %   from sample sd
          
          if ~exist('ev','var') || isempty(ev)
            ey=std(v_all_tmp,[],2)/sqrt(nv);
          else
            ey=ev/sqrt(nv);
          end
            
        else
            
            vmean=mean(v_all_tmp,1);
           v_all_tmp=v_all_tmp./repmat(vmean,[nt,1]);  % normalize 
           
           if ~exist('ev','var') || isempty(ev)
               
             y=mean(v_all_tmp,2);  
             ey=std(v_all_tmp,[],2)./sqrt(nv);
           
           else
               
           w= abs(repmat(vmean,[nt,1]));  % weight for weighted average
           
          y=sum(v_all_tmp.*w,2)./sum(w,2);
          
             ey=ev*mean(1./abs(vmean),2)/sqrt(nv);
           end
          
          
        end
    %    v_all_tmp=reshape(v_all_tmp,[size(b,ndims(b)),max(m(:))]);
     
        
        if j==1
            
        disp(tlt{i});
        fprintf(fid,'%s\n',tlt{i});
            fprintf('Mean velocity = %f (%f) cm/s\n',mean(v_all_tmp(:)),std(mean(v_all_tmp,1)));        
       
            fprintf(fid,'Mean velocity = %f (%f) cm/s\n',mean(v_all_tmp(:)),std(mean(v_all_tmp,1)));    
            vmean_all(i)=mean(v_all_tmp(:));
            
            %ind_max=floor(nt/2):floor(nt/2)+1;   
            if ~exist('ind_max','var')
                if mod(nt,2)==1
                    ind_max=[(nt-1)/2,(nt+1)/2,(nt+1)/2+1];
                else
                    ind_max=[nt/2,nt/2+1];
                end
                ind_max=2:nt-1;
            end
              
            if ~exist('ind_min','var')
                ind_min=[1,nt];
            end
            
          
            
            
          [pul,epul]=  pulsatility(v_all_tmp,ind_min,ind_max);  
          
          fprintf('Pulsatility = %f (%f)\n',pul,epul);
          fprintf(fid,'Pulsatility = %f (%f)\n',pul,epul);

          pulsa_all(i)=pul;
        end
      
        %   subplot(1,2,1);errorbar(1:size(b,ndims(b)),mean(v_all_tmp,2),std(v_all_tmp,[],2)/sqrt(max(m(:))));%
      
        if same_panel
           
          
            clr='rbg';
          hplot(end+1)=  errorbar((1:length(y))/nt,y,ey.*ones(nt,1),clr(i));%
          
        %    yl=[10,0];
            if i==1
                yl(1)=min(y-ey)-0.1;
                yl(2)=max(y+ey)+0.1;       
            else
                
                if exist('yl','var')
                yl(1)=min([yl(1),min(y-ey)-0.1]);
                yl(2)=max([yl(2),max(y+ey)+0.1]);
                else
                    yl(1)=min(y-ey)-0.1;
                    yl(2)=max(y+ey)+0.1;
                end
            end
            hold on;
            
            if i==2
                    ylim(yl);
            end
            if i==1
             plot((nshift+1)*[1,1]/nt,[yl(1),yl(1)+diff(yl)*0.1],'r-');
            end
        else
          
            errorbar((1:length(y))/nt,y,ey*ones(1,nt),'b');%
            
            %% plot pulse ox peak location
            
            
            ylim([min(y-ey)-0.1,max(y+ey)+0.1]);
            
            
            hold on;
            yl=ylim;
            
            plot((nshift+1)*[1,1]/nt,[yl(1),yl(1)+diff(yl)*0.1],'r-');
            
            title(tlt{i});
        end
        
        
        xlim([0.8,nt+0.2]/nt);
        
        set(gca,'TickLength',[0.02,0.02]);
        set(gca,'FontSize',12);
        if j==2
            ylabel('Normalized Velocity');
            xlabel('Cardiac Phase');
        else
            ylabel('Velocity (cm/s)');
            xlabel('Cardiac Phase');
           
        end
        
    end
   if j==1
       title(sprintf('v_{mean}=%3.2f; %3.2f cm/s',vmean_all(1),vmean_all(2)));
       legend(hplot,lgn);
   else
       title(sprintf('pulsa=%3.2f; %3.2f',pulsa_all(1),pulsa_all(2)));
   end
end

        fclose(fid);

