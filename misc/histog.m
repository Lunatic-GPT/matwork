function histog(dfiles,mask,xmin,xmax,nbin)
% histog(dfiles,mask,xmin,xmax,nbin)
% plot histogram of data in dfiles.
 
if ~iscell(dfiles)
    temp = dfiles;
    dfiles = cell(1);
    dfiles{1} = temp;
end

if ~exist('xmin','var') 
    xmin = 0;
end

if ~exist('nbin','var')
    nbin = 100;
end

nf = length(dfiles);
h = cell(1,nf);
figure; hold on;
symb = {'k-','b-','m-','dk','c*'};
ax_max = xmin;
for i=1:nf
    if ~exist('xmax','var')
    cmd = sprintf('3dhistog -mask ''%s'' -nbin %d -min %f ''%s'' > temp%d.1D', ...
                  mask,nbin,xmin,dfiles{i},i);
    else
    cmd = sprintf('3dhistog -mask ''%s'' -nbin %d -min %f -max %f ''%s'' > temp%d.1D', ...
                  mask,nbin,xmin,xmax,dfiles{i},i);
    end
    unix(cmd);
    h{i} = textread(['temp' num2str(i) '.1D'],'','commentstyle','shell');
    plot(h{i}(:,1),h{i}(:,2),symb{i});
    
    if ax_max < max(h{i}(:,1))
        ax_max = max(h{i}(:,1));
    end
    
   dfiles{i} = strrep(dfiles{i},'_','\_');    
end

mask = strrep(mask,'_','\_');

title(mask);
xlim([xmin,ax_max]);


legend(dfiles);