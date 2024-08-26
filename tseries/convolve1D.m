function convolve1D(hrfs,stim_pats,outf)
% convolve1D(hrfs,stim_pats,outf)
      nstim = length(hrfs);
      for i=1:nstim
         hrf = load(hrfs{i});
         pat = load(stim_pats{i});
         if i==1
             hrf_len = length(hrf);
             pat_len = length(pat);
             ts = zeros(hrf_len+pat_len-1,1);
         end
         
         hrf = reshape(hrf,hrf_len,1);
         pat = reshape(pat,pat_len,1);
         
         ts_conv = conv(hrf,pat);
         ts = ts+ts_conv;
         
      end
      plot(ts);
      save(outf,'ts','-ASCII');   
        
        