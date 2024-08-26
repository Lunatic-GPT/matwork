function tsroi(p)
% tsrois(xdata,ts_fnames,mask[,ts_scale,barplot,out_fname])
% ts_fnames: a cell array of dataset names.   
% ts_scale: an array to scale each dataset.
tic;

xdata = get(p,'x values');

ts_fnames = get(p,'data files');
temp = textscan(ts_fnames,'%s','Delimiter',',');
ts_fnames = temp{1};

mask = get(p,'roi mask');
temp = textscan(mask,'%s','Delimiter',';');
mask = maskcalc(temp{1}{:});

out_fname = get(p,'save file name');

if ~exist('barplot','var')
    barplot = 0;
end

nf = length(ts_fnames);

if  ~exist('ts_scale','var')|| isempty(ts_scale) 
    ts_scale =ones(1,nf);
end

fsize = 14;
for fi =1:nf
 
    cmd = sprintf('3dmaskave -mask ''%s'' ''%s'' > maskave_temp.1D',mask,ts_fnames{fi});
    unix(cmd);
    [ts_tmp,nvox_temp,temp_str] = textread('maskave_temp.1D','%f [%d %s');
    if ~strcmp(temp_str,'voxels]')
        error('3dmaskave output file error');
    end
    ts(:,fi) = ts_tmp(:)*ts_scale(fi);
end
ts_len = size(ts,1);
nvox = nvox_temp(1);
  figure;
  if barplot
      h = bar(ts);
      set(h,'BarWidth',1);
      colormap([0.6,0.6,0.6;0,0,1;0.6,0,0]);
  else
       hold on;
      % symb = {'k+','bs','r*','g>'};
       % symb = {'ks','b*','r+','g>','m','c'};
        symb = {'k','b','r','g','m','c'};    
        %symb = {'k-','b*','rs','k-','k-'};
       for i=1:size(ts,2)
        plot(xdata,ts(:,i),symb{i});
       end
  end
  
  legend('toggle');
  down = min(ts(:));
  up = max(ts(:));
  nvox_str =sprintf('n=%d',nvox);
 text(ts_len/2,down+(up-down)/4,nvox_str,'FontSize',fsize);
 mask = strrep(mask,'_','\_');
 title(mask);
 set(gca,'box','off');
 
 %xlim([0,15]);
 %ylim([0,up*1.1]);
if exist('out_fname','var')
    params = p;
  save(out_fname,'ts','mask','ts_fnames','nvox','params');
end

       
disp([mfilename ' finish in ', num2str(toc), ' s']);
     
    


      
         
         
         
         


