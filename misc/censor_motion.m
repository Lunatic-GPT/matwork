function censor_motion(t1,t2,TR)
% censor_motion(t1,t2,TR)
% t1: the translation velocity threshold (mm/s)
% t2: the rotation velocity threshold (degree/s)
if ~exist('t1','var')
    t1 = 0.05/2;
end
if ~exist('t2','var')
    t2 = 0.05/2;
end
if ~exist('TR','var')
    TR = 2;
end

a =dir('motion_s*.1D');
fid = fopen('censor_points.1D','w');
for i=1:length(a)
    
    d = load(a(i).name);
    
    s = abs(diff(d,1,1)/TR);
    
    c = (s(:,1)>t1 | s(:,2)>t1 | s(:,3)>t1 | ...
         s(:,4)>t2 | s(:,5)>t2 | s(:,6)>t2);
     
     in = find(c);
     in = [in,in+1];
     in = unique(in);
     if ~isempty(in)
         fprintf(fid,'%d,',in-1);
         fprintf(fid,'\n');
     else
         fprintf(fid,'-1\n');
     end
end