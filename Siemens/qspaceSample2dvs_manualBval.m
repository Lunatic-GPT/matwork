function qspaceSample2dvs_manualBval(fname,bzero,bvals,repetitions,b_order,invert_vec)
% bzero: number of bzero.
% b_order: one based

if ~exist('invert_vec','var')
    invert_vec=0;
end

a=load(fname);

nshell=length(unique(a(:,1)));

if nshell~=length(bvals)
   error('Number of elements in bvals is wrong');
end


if ~exist('bzero','var')
    bzero=round(size(a,1)/nshell/6);
end

  if min(a(:,1))==0
      a(:,1)=a(:,1)+1;
  end
v=a(:,2:4);

maxb=max(bvals);
for i=1:nshell
      v(i==a(:,1),:)=v(i==a(:,1),:)*sqrt(bvals(i)/maxb);
end

v2=v*0;
start=1;
for i=1:nshell
   
    n=sum(b_order(i)==a(:,1));
    v2(start:start+n-1,:)=v(b_order(i)==a(:,1),:);
    
    start=start+n;
    
    
end





step=ceil(size(v,1)/bzero);

for i=1:bzero
    
    if i==1
       v3=cat(1,[0,0,0],v2((i-1)*step+1:i*step,:)); 
    elseif i<bzero
     v3=cat(1,v3,[0,0,0],v2((i-1)*step+1:i*step,:));
    else
    v3=cat(1,v3,[0,0,0],v2((i-1)*step+1:end,:));
    end    
    
end

v3=repmat(v3,[repetitions,1]);
    
prefix=strtok(fname,'.');


prefix=[prefix,sprintf('_%d',bvals),'_rep',num2str(repetitions)];

if invert_vec
    v3=-v3;
    prefix=[prefix,'_inv'];
end


write_DiffusionVectors(v3,prefix);



