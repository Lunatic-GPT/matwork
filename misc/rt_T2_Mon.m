function rt_T2_Mon(fid_prefix,exp_dir)
%rt_T2_Mon(fid_dir,exp_dir)

m=get_roi_mask();
%m=zeros(64,64,4);
%m(32,32,2)=1;
seqfil='obl_epi8_MT';
if ~exist('fid_prefix','var')
fid_prefix = 'acqfil';
end
%exp_dir = '/home/xiaopeng/vesta/vnmrsys/exp1/';
if ~exist('exp_dir','var')
exp_dir='/data/users/xiaopeng/vnmrsys/exp1';
end
%exp_dir=pwd;

yl = [0,18.4]*1e4;
fid_dir=fullfile(exp_dir,fid_prefix);
hg= figure;
datenum_old =0;
n_proc=0;
   te=readPar(fid_dir,'te');
    te=te(2:end);
    sf=readPar(fid_dir,'seqfil');
    
    if ~strcmp(sf(2:end-1),seqfil)
        error('seqfil mismatch');
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
      plot(te(1:length(ts)),log10(ts),'o');
      xlabel('Time (s)');
      
    else 
      x=get(hl,'XData');
      y=get(hl,'YData');
      y=[y,log10(ts)];
      set(hl,'YData',y);
      set(hl,'XData',te(1:length(ts)+length(x)));
      
    end
    xlim([0,0.07]);
%    ylim(yl);

    pause(5);
end

    