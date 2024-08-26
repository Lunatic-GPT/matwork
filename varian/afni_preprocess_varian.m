function afni_preprocess_varian(flist)
% afni_preprocess_varian(flist)
% data preprocessing in afni.
% flist: afni data files, use subbrik selection to exclude dummy scans.


flist = str2cell(flist);
nf = length(flist);

flist2=cell(1,nf);
for i=1:nf
   if ~any(flist{i}=='+')
       flist{i}=[flist{i},'+orig'];
   end
    
  prefix_out = ts_dtr(flist{i});
    
    prefix_out2 = mean_signal(flist{i});
    
    percent_change([prefix_out,'+orig'],[prefix_out2,'+orig']);
    
    flist2{i} = [prefix_out2,'+orig'];
end

 brainmask(flist2,0.55);

 
function prefix_out=mean_signal(fname)

    prefix = strtok(fname,'+');
    cmd =  sprintf('3dTstat -prefix %s_mean %s',prefix,fname);
    prefix_out = [prefix,'_mean'];
    unix(cmd); 
    
    
function brainmask(flist,mfrac)

nf = length(flist);
for i=1:nf 
    
    cmd = sprintf('3dClipLevel -mfrac %f %s',mfrac,flist{i});
    [s,clip_tmp] = unix(cmd);
    clip(i)=str2double(clip_tmp);
       
end

cmd='3dcalc ';
expr = '';
for i=1:nf
    cmd = sprintf('%s -%c %s',cmd,'a'+i-1,flist{i});
    expr = sprintf('%sstep(%c-%f)*',expr,'a'+i-1,clip(i));
end

unix([cmd,' -expr ''',expr(1:end-1),''' -prefix brainmask']);
    

function percent_change(f_ts,f_mean)
    
    prefix = strtok(f_ts,'+');
    cmd = sprintf('3dcalc -a %s -b %s -expr ''a/b'' -prefix %s_norm',f_ts,f_mean,prefix);
    unix(cmd);
   
    


function prefix_out = ts_dtr(f_in)

prefix=strtok(f_in,'+');
cmd = sprintf('3dDetrend -prefix %s_dtr -polort 2 %s',prefix,f_in);
unix(cmd);
prefix_out = [prefix,'_dtr'];
   
function prefix_out = low_pass(f_in,freq_thr)   

prefix=strtok(f_in,'+');

   cmd = sprintf('3dFourier -prefix %s_lp -lowpass %f %s',prefix,f_in,freq_thr);
   unix(cmd);
   
prefix_out = [prefix,'_lp'];

