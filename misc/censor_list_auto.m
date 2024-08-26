function censor_list_auto(fpat,nthr,mask)
% censor_list_auto(fpat,nthr,mask)

if ~exist('mask','var')
    mask = 'brainmask+orig';
end

if ~exist('nthr','var')
   nthr = 50;
end

fpat = strtok(fpat,'.');

a =dir([fpat,'.HEAD']);
fid = fopen('censor_points.1D','w');
for i=1:length(a)

    cmd = sprintf('3dToutcount -mask %s %s > 3dToutcount_out.txt',mask,a(i).name);
    unix(cmd);
    nout = load('3dToutcount_out.txt');
    in = find(nout>=nthr);   
     if ~isempty(in)
         fprintf(fid,'%d,',in-1);
         fprintf(fid,'\n');
     else
         fprintf(fid,'-1\n');
     end
end