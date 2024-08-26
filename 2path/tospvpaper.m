function dname=tospvpaper
d=mfilename('fullpath');
d=fileparts(d);
cd(fullfile(d,'..','papers_abstracts','partialVolume_SWI')); 
dname=fullfile(d,'..','papers_abstracts','partialVolume_SWI');