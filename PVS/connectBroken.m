
function connectBroken(prefix)
%%
fprintf('not tested yet');

disp('');
a=load(prefix);
%subind=a.subind;

voxsize=a.voxsize;
i=1;
pvs_combined=0*a.c;
while i<length(a.ind)
   
    combined=false;
    for j=i+1:npvs
        x=a.subind{i}.*voxsize;
        y=a.subind{j}.*voxsize;
        [res,pair]=shouldCombine(x(1:2,:),y(1:2,:));       
       if  res
          if ~any(pair~=[-1,1])
               [a,pvs_combined]=do_combine(a,i,j,pair,pvs_combined); 
              combined=true;
              break;
          else
              if a.nvox_path(i)>3
                 error('PVS %d - %d: combine test error', i, j);
              end
          end
          
       end
          
       
        [res,pair]=shouldCombine(x(1:2,:),y(end-1:end,:));       
       if  res
          if ~any(pair~=[-1,-1])
               [a,pvs_combined]=do_combine(a,i,j,pair,pvs_combined); 
              combined=true;
              break;
          else
              if nvox_path(i)>3
                 error('PVS %d - %d: combine test error', i, j);
              end
          end
          combine_;
       end
       
         [res,pair]=shouldCombine(x(end-1:end,:),y(1:2,:));       
       if  res
          if ~any(pair~=[1,1])
               [a,pvs_combined]=do_combine(a,i,j,pair,pvs_combined); 
              combined=true;
              break;
          else
              if nvox_path(i)>3
                 error('PVS %d - %d: combine test error', i, j);
              end
          end
          combine;
       end
       
        [res,pair]=shouldCombine(x(end-1:end,:),y(end-1:end,:));       
       if  res
          if ~any(pair~=[1,-1])
                  [a,pvs_combined]=do_combine(a,i,j,pair,pvs_combined);   
                  combined=true;
                  break;
          else
              if nvox_path(i)>3
                 error('PVS %d - %d: combine test error', i, j);
              end
          end
        
       end

    
%            
%            ...
%            || shouldCombine(x(end-1:end,:),x(1:2,:)) ...
%            || shouldCombine(x(1:2,:),x(end-1:end,:)) ...
%            || shouldCombine(x(end-1:end,:),x(1:2,:))
%       
%            pair(end+1,:)=[i,j];
%               
      
    end
    
    if ~combined %if not combined go to next PVS, otherwise check again.
      i=i+1;
    end
  
end
 
function [a,pvs_combined]=do_combine(a,i,j,pair,pvs_combined)

fprintf('do_combine not finished');
if pair(1)==1 && pair(2)==1
  a.nvox_path(end+1)=sum(a.nvox_path([i,j]));
  a.ind{end+1}=[a.ind{i}(1:n1),a.ind{j}(1:n2),a.ind{i}(n1+1:end),a.ind{j}(n2+1:end)];
  a.subind{end+1}=cat(1,a.subind{i}(1:n1,:),a.subind{j}(1:n2,:),a.subind{i}(n1+1,:),a.subind{j}(n2+1,:));
  a.c(a.ind{end+1})=length(a.nvox_path)-2;
  
elseif pair(1)==1 && pair(2)==-1
    
elseif pair(1)==-1 && pair(2)==1
    
elseif pair(1)==-1 && pair(2)==-1
    
end  

a.nvox_path([i,j])=[];
a.ind([i,j])=[];
a.subind([i,j],:)=[];
val=pvs_combined(a.ind{i}(1));

if val>0  %if PVS i already registered in the pvs_combined;
   pvs_combined(a.ind{i}(1))=value;  
else
   pvs_combined(a.ind{i}(1))=max(pvs_combined(:))+1; 
end




function [res,pair]=shouldCombine(x1,x2)

% if [1,1]: the combined path should be [x1,x2] 
% if [1,-1]: the combined path should be [x1, flip(x2)];
% if [-1,1]: the combined path should be [flip(x1),x2]
% if [-1,-1]: the combined path should be [flip(x1),flip(x2)]
dx1=x1(1,:)-x1(2,:);
dx1=dx1/sos(dx1);
dx2=x2(1,:)-x2(2,:);
dx2=dx2/sos(dx2);

cs=sum(dx1.*dx2);
flipped=false;
if cs<0
    dx2=-dx2;
    x2=flip(x2,1);
    flipped=true;
end

thr=cos(10*pi/180);
thr_dist=0.83;
res=false;
x1_to_x2=true;
if (cs>thr)
    [d,t]=minDistance(x1(2,:),dx1,x2(1,:));
    
    if t>0
        if d<thr_dist
          res=true;
        end
    else
        [d,t]=minDistance(x2(2,:),dx2,x1(1,:));
        if t>0 && res<thr_dist
            res=true;
            x1_to_x2=false;
        end     
    end
  
end


pair=[];
if res
  if x1_to_x2
      if flipped
          pair=[1,-1];
      else
          pair=[1,1];
      end
  else
      if flipped
          pair=[-1,1];
      else
          pair=[-1,-1];
      end
      
  end
end

