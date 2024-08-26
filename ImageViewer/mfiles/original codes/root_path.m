function root_path
% root_path defines the main directory for the different .m-files that are
% used by BLAT

root = uigetdir(matlabroot,'Choose the BLAT Root Directory');

save root_directory root

f1=fullfile(root,'temporary_saved_variables');
if (exist(f1) == 0)
   mkdir (f1);
end