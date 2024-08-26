function gating_resize_file(fpattern,dummy,dur)
%gating_resize_file(fpattern,dummy,dur)
% dummy,dur: unit second.
% dur: duration of the fMRI scan in addition to dummy scan.
nd = dummy*40;
n = dur*40;

dir_str = dir(fpattern);
for i=1:length(dir_str)
    a = load(dir_str(i).name);
    b = a(nd+1:nd+n);
    prefix = strtok(dir_str(i).name,'.');
    save_mat_int(b,[prefix,'_ndum.1D']);
    
    fprintf('Discard the last %3.2f second in %s\n',(length(a)-nd-n)/40,dir_str(i).name);
end

