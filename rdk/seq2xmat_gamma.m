function seq2xmat_gamma(sti_file,nTR_pre,nTR_off1,nTR_on,nTR_off2,TR,xmat_file)
% seq2xmat_gamma(sti_file,nTR_pre,nTR_off1,nTR_on,nTR_off2,TR,xmat_file)
% sti_file: a 2D matrix of stimulus sequences.  Each row corresponds to a
% different scan.  
% nTR_pre: number of baseline TRs at the beginning of each scan
%nTR_off1: baseline TR before stimulus in each trial.
%nTR_on: duration for stimulus on
%nTR_off2: baseline TR after stimulus in each trial.
% TR in second.

if ~exist('TR','var')
    TR = 2;
end

if ~exist('power','var')
    power = 8.6;
end

if ~exist('time_delay','var')
    time_delay = 0;
end

if ~exist('width','var')
    width = 0.547;
end

if isa(sti_file,'char')
  sti_seq = load(sti_file);
else
  sti_seq = sti_file;
end

npat = max(sti_seq(:));
nTR_scan = nTR_pre+(nTR_off1+nTR_on+nTR_off2)*size(sti_seq,2);
xmat = zeros(nTR_scan*size(sti_seq,1),npat);    

%% generate gamma function in data
           data = [];
           t = 0;
           peak = (power*width)^power*exp(-power); % peak value of the gamma function.
           while (t<time_delay+power*width || data(end)>0.01*peak)
             l = t-time_delay;
             if l>0
             data(end+1) = (l^power)*exp(-l/width);
             else
             data(end+1)=0;
             end
             t = t+TR;
           end
           
for j=1:size(sti_seq,1) % loop through scans
  
arr = zeros(nTR_scan,npat);

for i = 1:size(sti_seq,2)  %loop over the different trial types

    ind1 = nTR_pre+(nTR_on+nTR_off1+nTR_off2)*(i-1)+1;
    ind2 = ind1+nTR_on-1;
    arr(ind1:ind2,sti_seq(j,i)) = 1;
end
     arr_gamma = conv2(data(:),1,arr);
     xmat(nTR_scan*(j-1)+1:nTR_scan*j,:) = arr_gamma(1:nTR_scan,:);   

end
      


fid=fopen(xmat_file,'w');

for i=1:size(xmat,1)
    fprintf(fid,'%4.2f ',xmat(i,:));
    fprintf(fid,'\n');
end

fclose(fid);