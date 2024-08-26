fid=fopen('c:/temp/passtemp.txt','r');
passwd=textscan(fid,'%s');
passwd=passwd{1}{1};
fclose(fid);

sid={'HVB267','HVB268','HVB269','HVB270','HVB271','LTL0148','LTL0156', 'LFT0157','LTL0158', 'RTL0242', 'RTL0159', 'RTL0255'};

for i=2:length(sid)
    disp(sid{i});
cd(sid{i});    
src=sprintf('/netscr/zongx/CEST_multiPool/%s/run_dofit_scan_1_2_fixoffset_fixpar_1_2_9_11_17_18_pixbypix_combined_T1w5.0_Consistent_TI_B1map.mat',sid{i});
 %src=sprintf('/netscr/zongx/CEST_multiPool/%s/subjobs_2c.sh',sid{i});
    
status=scpfrommatlab('zongx','killdevil.unc.edu','localhost',passwd,src,'dofit_Results_temp.mat');
if status==0
    disp('OK');
else
    disp('Transfer Failed');
end
cd ..;
end


