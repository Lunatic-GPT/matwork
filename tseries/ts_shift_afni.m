function data = ts_shift_afni(data,ind_pts)
%data = ts_shift_afni(data,ind_pts)
% use 3dcalc in afni to shift the data so that the mean of ind_pts are
% zero. 
% ind_pts is a string, zero based and has the same format as in afni such
% as '[0,1,38..42]'

prefix = strtok(data,'+');
cmd = sprintf('3dTstat -mean -prefix %s_bl %s''%s''',prefix,data,ind_pts);   
unix(cmd);

cmd = sprintf('3dcalc -prefix %s_shft -a %s -b %s_bl+orig -expr ''(a-b)''',prefix,data,prefix);
unix(cmd);
        


