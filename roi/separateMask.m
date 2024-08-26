function separateMask(dname,cs,direction,out_prefix)
% separateMask(dname,cs,direction,out_prefix)
% dname: mask name
% cs: minimum cluster size
% direction: 'lr','ap', or 'is';
% out_prefix: output mask with two values.  1 corresponds to smaller index,
% 2 larger indices.  Whether 1 is left, anterior, or inferior will depend
% on the orientation code.  
maskcalc(dname,'a',cs,'separateMask_temp');
if isempty(findstr(dname,'+tlrc'))
    [d,info] = BrikLoad('separateMask_temp+orig');
else
    [d,info] = BrikLoad('separateMask_temp+tlrc');
end
nc = max(d(:));
d1=zeros(size(d));
%orient = info.Orientation;

for i=1:nc
    clus = find(d==i);
    pos = zeros(1,length(clus));
    for j=1:length(clus)
        [i1,i2,i3] = ind2sub(size(d),clus(j));
        switch direction
            case 'lr'
              pos(j) = i1;
            case 'ap'
                pos(j) = i2;
            case 'is'
                pos(j) = i3;
            otherwise
                error('wrong value of direction');
        end
    end
 
    [pos_temp,ind] = sort(pos);
    
    for j=1:floor(length(pos)/2)
        [i1,i2,i3] = ind2sub(size(d),clus(ind(j)));
        d1(i1,i2,i3) =1;
    end
        
end

d2=2*(d>0 & (d1 == 0));
data = d1+d2;
Opt.OverWrite = 'y';
if ~isempty(strfind(dname,'+tlrc'))
  Opt.View = '+tlrc';
end
 
Opt.Prefix = out_prefix;
WriteBrik(data,info,Opt);

    
    
