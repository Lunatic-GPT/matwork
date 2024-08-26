function separateMasktlrc(dname,tlrc,out_prefix)
% separateMasktlrc(dname,tlrc,out_prefix)
% dname: mask name
% tlrc: name of tlrc dataset
% out_prefix: output mask with two values.  1 corresponds to smaller index,
a1 = tlrc2orig(0,0,0,tlrc);
a2 = tlrc2orig(0,1,1,tlrc);
a3 = tlrc2orig(0,1,-1,tlrc);

m = cross(a1-a2,a1-a3);
m = m/sqrt(sum(m.*m));


[d,info] = BrikLoad(dname);
lr = zeros(size(d,1),size(d,2),size(d,3));

for i=1:size(d,1)
    for j=1:size(d,2)
        for k=1:size(d,3)         
            x = info.ORIGIN+info.DELTA.*[i-1,j-1,k-1];
            if (x-a2')*m > 0
                lr(i,j,k)=1;
            else 
                lr(i,j,k) = 2;
            end
        end
    end
end

lr(d==0) = 0;
Opt.Prefix = out_prefix;
history = sprintf('separateMasktlrc(%s,%s,%s)',dname,tlrc,out_prefix);
WriteBrikEZ(lr,info,history,out_prefix);

    
    
