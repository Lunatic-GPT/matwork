function clist = tsAlign_censor(cfile,sti_seq,frames,nTR_trial)
%clist = tsAlign_censor(cfile,sti_seq,frames)
% cfile: file containing censored time points
% sti_seq: stimuus sequence file
% frames: time frames in each scan
% allow different durations of different trials.  
a= read_censor_points(cfile);

seq = load(sti_seq);
seq = seq';
seq = seq(:);
ns = length(a);
censor = zeros(1,ns*length(frames));

for i=1:ns
    for j=1:length(a{i}) 
      ind = find(a{i}(j)==frames);
      if ~isempty(ind)
        censor(ind+(i-1)*length(frames)) = 1;
      end
    end
end

npat = max(seq(:));
clist = cell(1,npat);
if ~exist('nTR_trial','var')
  nTR_trial = ns*length(frames)/length(seq)*ones(1,npat);
end

for i=1:npat
    ntrl = length(find(seq == i));
    clist{i} = zeros(1,nTR_trial(i)*ntrl);
end

ind = find(censor>0);

Tfin = cumsum(nTR_trial(seq)); % finish time of the trials

for i=1:length(ind)
    
    itrl = find(ind(i)>Tfin(1:end-1) & ind(i) <=Tfin(2:end));
    
    if isempty(itrl)
        itrl = 1;
        os = ind(i);
    else
        os = ind(i) - Tfin(itrl);
        itrl = itrl+1;
    end
    
    ipat = seq(itrl);
    
    itrl_pat = length(find(seq(1:itrl)==ipat));
    clist{ipat}((itrl_pat-1)*nTR_trial(ipat)+os) = 1;

end
    