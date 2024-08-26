function tsAlign_gui

p = parameter('sort time series from different trials');

p = add(p,'string','time series','');
p = add(p,'button','Browse','browse_callback(params,''time series'')');
p = add(p,'string','subbrik range','0..100'); 
p = add(p,'filename','stimulus sequence file','');

p = add(p,'int','trial duration (TR)',17);

p = add(p,'int','baseline duration (TR)',0);  
% the last few points will be used as baseline to adjust the time series trial by trial.
% to skip the adjustment, set the baseline value to 0 or empty.

p = add(p,'string','prefix for output files','ts_[filter name]');
p = add(p,'button','Go','tsAlign_Sepfile(params)');

p = add(p,'button','generate sequence file','eprime_seq_gui');
p = add(p,'button','Quit','close');

parametergui(p);