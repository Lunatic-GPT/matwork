function Layers_Between_Boundaries(fname,nl)
%  Layers_Between_Boundaries(fname,nl)
% fname: name of the make containing the boundary definitions with values
% equals to 1 and 2, respectively.
% nl:  number of layers to generate.
% only works for single slice data.
[a,info] = BrikLoad(fname);

mask = zeros(size(a));
for islc = 1:size(a,3)
    
[x1,y1]=find(a(:,:,islc)==1);
[x2,y2]=find(a(:,:,islc)==2);
    
n1=length(x1);
n2=length(x2);

if n1<1 || n2<1
    continue;
end

c1r=order([x1,y1]);
c2r=order([x2,y2]);

% make 1 longer than 2
if n1<n2
    tmp=n1;
    n1=n2;
    n2=tmp;
    tmp=c1r;
    c1r=c2r;
    c2r=tmp;
end

grid=zeros(2,n1,nl+1);
% generate the grid x
for i=1:n1
    
    i2=1+round((n2-1)*(i-1)/(n1-1));
    grid(1,i,:) = linspace(c1r(i,1),c2r(i2,1),nl+1);
    grid(2,i,:) = linspace(c1r(i,2),c2r(i2,2),nl+1);
    
end

% loop through the grids and find the voxel that
% fall within the layers.

x_max = max([x1;x2]);
x_min = min([x1;x2]);
y_max = max([y1;y2]);
y_min = min([y1;y2]);

area = cell(x_max-x_min+1,y_max-y_min+1);
grid_ind = cell(x_max-x_min+1,y_max-y_min+1);

for i=1:nl
    for j=1:n1-1
        
        [area_tmp,ind_tmp,nvox] = intersect_voxels(grid(:,j:j+1,i:i+1));
        for k=1:nvox
          u=ind_tmp(k,:)-[x_min,y_min]+1;
          grid_ind{u(1),u(2)}(end+1,:) = [j,i];
          area{u(1),u(2)}(end+1) = area_tmp(k);
        end
       
    end
end

for i=1:x_max-x_min+1
    for j=1:y_max-y_min+1
      if isempty(area{i,j})
          continue;
      end
      [tmp,ind]=max(area{i,j});
      mask(i+x_min-1,j+y_min-1,islc) = grid_ind{i,j}(ind,2);
    end
end

end

prefix = strtok(fname,'+');
prefix = [prefix,'_layers'];
history = sprintf('Layers_Between_Boundaries(%s,number of layers = %d)',fname,nl);
WriteBrikEZ(mask,info,history,prefix);


function [area,ind,n]= intersect_voxels(square)
% given the coordinates of 4 vertices of a square, find the voxels that intersect with it.
% Voxel [i,j] has vertices of [i+-0.5,j+-0.5].

xs=square(1,:,:);
xs=xs(:);
ys=square(2,:,:);
ys=ys(:);

xmin=round(min(xs));
xmax=round(max(xs));

ymin=round(min(ys));
ymax=round(max(ys));
 n=0;
 ind = [];
 area = [];
for x=xmin:xmax
    for y=ymin:ymax
        
      xv=[x-0.5,x-0.5,x+0.5,x+0.5];
      yv=[y-0.5,y+0.5,y-0.5,y+0.5];
      
      [Dx,Dy] = polybool('intersection',xv,yv,xs,ys);
      
      if ~isempty(Dx)
          area(end+1) = polyarea(Dx,Dy);
          n=n+1;
          ind(end+1,:) =[x,y];
      end
        
    end
end


function x2=order(x)

% order the points in x such that they are in the same order as in a curve.
% each point should have two or less nearest neighbors and two or less
% next-nearest neighbors.
np = size(x,1);
nn = zeros(1,np);
nnn=zeros(1,np);
indn=cell(1,np);
indnn=cell(1,np);
ind_terminal = [];
for i=1:np
    
    [nn(i),indn{i}] = nearest_neighbors(x(i,:),x);
    [nnn(i),indnn{i}] = next_nearest_neighbors(x(i,:),x);   
    if nn(i)>2 
     error('Points have more than two nearest neighbors');
    end
    if nnn(i)>2 
     error('Points have more than two nearest neighbors');
    end
    if nn(i)+nnn(i)<1
     error('Points have no nearest and next-nearest neighbors');
    end

    if is_terminal(x(i,:),x(indn{i},:),x(indnn{i},:))
      ind_terminal(end+1) = i;
    end
end



if length(ind_terminal)~=2
    error('Number of terminals %d',length(ind_terminal));
end

ind_ordered(1) = ind_terminal(1);

while length(ind_ordered)<np
    
    
    indn_tmp = setdiff(indn{ind_ordered(end)},ind_ordered);
    indnn_tmp = setdiff(indnn{ind_ordered(end)},ind_ordered);
    
    
    if length(indn_tmp)==1
        ind_ordered(end+1) = indn_tmp;
    elseif length(indnn_tmp)==1
        ind_ordered(end+1) = indnn_tmp;
    else
        error('Order error');
    end
end

x2 = x(ind_ordered,:);
    


%first find a terminal



function [n,ind] = nearest_neighbors(x,xlist)

n=0;
ind = [];
  for i=1:size(xlist,1)
      
      if x(1) == xlist(i,1)+1 && x(2) == xlist(i,2)
         n=n+1;
         ind(end+1) = i;
      elseif x(1) == xlist(i,1)-1 && x(2) == xlist(i,2)
          n=n+1;
         ind(end+1) = i;
      elseif x(1) == xlist(i,1) && x(2) == xlist(i,2)+1
          n=n+1;
         ind(end+1) = i;
      elseif x(1) == xlist(i,1) && x(2) == xlist(i,2)-1
         n=n+1;
         ind(end+1) = i;
      else
          
      end
  end

function [n,ind] = next_nearest_neighbors(x,xlist)

n=0;
ind = [];
  for i=1:size(xlist,1)
      
      if x(1) == xlist(i,1)+1 && x(2) == xlist(i,2)+1
         n=n+1;
         ind(end+1) = i;
      elseif x(1) == xlist(i,1)-1 && x(2) == xlist(i,2)-1
          n=n+1;
         ind(end+1) = i;
      elseif x(1) == xlist(i,1)-1 && x(2) == xlist(i,2)+1
          n=n+1;
         ind(end+1) = i;
      elseif x(1) == xlist(i,1)+1 && x(2) == xlist(i,2)-1
         n=n+1;
         ind(end+1) = i;
      else
          
      end
  end
  
 function isterminal = is_terminal(x,xn,xnn)
     
     isterminal = true;
     if size(xn,1)==2 || size(xnn,1)==2
         isterminal = false;
     end
     
     if size(xn,1)==1 && size(xnn,1)==1
         if (x(1)-xn(1))*(x(1)-xnn(1))<0
             isterminal = false;
         end
         
         if (x(2)-xn(2))*(x(2)-xnn(2))<0
             isterminal = false;
         end
     end
     
         
  
  
  
      
      
      
      








