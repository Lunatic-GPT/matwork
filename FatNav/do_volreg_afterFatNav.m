function do_volreg_afterFatNav(mid,s,center)

% center in the [ro,pe,par] in the dicom convention.

if ~exist('center','var')
    center = [0,0,0];
end

    wait_for_files([mid,'_FatNav'],s);

     mat2afni_FatNav([mid,'_FatNav'],s,center);
     
cmd=sprintf('3dvolreg -dfile Motion_%s.1D -base 0 -prefix %s_FatNav_volreg %s_FatNav+orig',mid,mid,mid);
disp(cmd);
unix(cmd);
sid=filename(pwd);
mkdir(sprintf('~/MotionFiles/%s',sid));

cmd=sprintf('cp Motion_%s.1D ~/MotionFiles/%s',mid,sid);  %the file in netscr cannot be found on longleaf jobs created with sbatch although can be found inside matlab
disp(cmd);
unix(cmd);


function wait_for_files(prefix,ind1,nsecond)
found =false;

if ~exist('nsecond','var')
    nsecond = 60;
end

while ~found
    pause(nsecond);
    found=true;
for i=1:length(ind1)
  
    name=sprintf('%s_%d_*.mat',prefix,ind1(i));
    dir_str=dir(name);
    
    if isempty(dir_str) || dir_str(1).bytes<1000
        found = false;
        break;
    end

end

end