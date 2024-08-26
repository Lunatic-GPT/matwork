function showsScans_dcm(maxScan)

if ~exist('maxScan','var')
    maxScan=150;
end


for i=1:maxScan
    
    if exist(num2str(i),'dir') 
    
        str=dir([num2str(i),'\*.dcm']);
     a=dicominfo(fullfile(num2str(i),str(1).name));
     
         fprintf('%d %s; ',i, a.ProtocolName);
     if isfield(a,'ImageComments')
     %  fprintf('%s;  ',a.ImageComments);
     end
     
     if isfield(a,'SeriesDescription')
         fprintf('%s;', a.SeriesDescription);
     end
     fprintf('\n');
    end
    
end
