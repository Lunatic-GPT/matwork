m=get_roi_mask();
%m=zeros(64,64,4);
%m(32,32,2)=1;
seqfil='obl_epi8_MT';
fid_prefix = 'acqfil';
%exp_dir = '/home/xiaopeng/vesta/vnmrsys/exp1/';
exp_dir='/data/users/xiaopeng/vnmrsys/exp7';
%exp_dir=pwd;

yl = [1.5,3.5]*1e4;
fid_dir=fullfile(exp_dir,fid_prefix);
hg= figure;
datenum_old =0;
n_proc=0;

    tr=readPar(fid_dir,'total_tr');
    tr=tr/1000;
    sf=readPar(fid_dir,'seqfil');
    if ~strcmp(sf(2:end-1),seqfil)
        error('seqfile mismatch');
    end
    
  while(1)
        
    dir_str=dir([fid_dir,'/fid']);
    if isempty(dir_str)
        continue;
    end
    if datenum_old~=dir_str.datenum
    datenum_old = dir_str.datenum;
    [epi,n_proc]=epiShaper_rt(fid_dir,n_proc);
    else 
        continue;
    end
    
    if isempty(epi)
        continue;
    end
    ts=mean_roi(epi,m);
    
    figure(hg);
    hl=findobj(hg,'Type','line');
    if isempty(hl)
      plot(tr*(1:length(ts)),ts);
      xlabel('Time (s)');
      if length(ts)<180
          xlim([0,180]);
      else
          xlim([length(ts)-180,length(ts)]);
      end
   
    else 
      x=get(hl,'XData');
      y=get(hl,'YData');
      y=[y,ts];
      set(hl,'YData',y);
      set(hl,'XData',tr*(1:length(ts)+length(x)));
      xl=xlim;
      if tr*(length(ts)+length(x))>xl(2)
        xl=xl+tr*length(ts);
        xlim(xl);
      end
        ind = round(xl/tr);
        if ind(2)>length(y)
            ind(2)=length(y);
        end
        if ind(1)<1
            ind(1)=1;
        end
    %    yl = [min(y(ind(1):ind(2))),max(y(ind(1):ind(2)))];  
      
    end
    ylim(yl);
    
    pause(5);
end

    