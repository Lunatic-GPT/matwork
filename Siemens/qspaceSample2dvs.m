function qspaceSample2dvs(fname,bzero,uniform_q)
% bzero: number of bzero.

if isa(fname,'char')
    a=load(fname);
else
    a=fname;
end


nshell=length(unique(a(:,1)));

if ~exist('bzero','var')
    bzero=round(size(a,1)/nshell/6);
end


v=a(:,2:4);

for i=1:nshell

    if min(a(:,1))==1
        if uniform_q
            v(i==a(:,1),:)=v(i==a(:,1),:)*(i/nshell);
        else
            v(i==a(:,1),:)=v(i==a(:,1),:)*sqrt(i/nshell);
        end
    else  % zero based
        if uniform_q
            v(i-1==a(:,1),:)=v(i-1==a(:,1),:)*(i/nshell);
        else
            v(i-1==a(:,1),:)=v(i-1==a(:,1),:)*sqrt(i/nshell);
        end
    end

end

step=ceil(size(v,1)/bzero);

for i=1:bzero
    
    if i==1
       v2=cat(1,[0,0,0],v((i-1)*step+1:i*step,:)); 
    elseif i<bzero
     v2=cat(1,v2,[0,0,0],v((i-1)*step+1:i*step,:));
    else
    v2=cat(1,v2,[0,0,0],v((i-1)*step+1:end,:));
    end    
    
end

if bzero==0
    v2=v;
end

if isa(fname,'char')    
prefix=strtok(fname,'.');
else
    
    prefix='qspaceSampe2dvs';
end

if uniform_q
    prefix=[prefix,'_uniformQ'];
else
    prefix=[prefix,'_uniformB'];
end

write_DiffusionVectors(v2,prefix);



