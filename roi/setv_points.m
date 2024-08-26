function res=setv_points(mask,points,val)
% res=setv_points(mask,points,val)
%points: n*3; 
% val: set points to val in mask;

res=mask;
for i=1:size(points,1)
    res(points(i,1),points(i,2),points(i,3))=val;
   
end
