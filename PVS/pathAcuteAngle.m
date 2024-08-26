
%[p1,p2,p3]=
pathAcuteAngle(prefix2);  % detect PVS that forms an acute angle along % its path; 
                          % do not implement for now; only 1 - 2 PVSs has acute angles 
                          % and hard to separate 1 PVS into 2 PVSs in such cases
                          

function pathAcuteAngle(prefix)
% file should contain a ind file (1*npath)
% path file (1*npath)
% c: non zero with values equal to the npath index.
% remove path with >90 degree turns, i.e. forming acute angles

global voxsize;
a=load([prefix,'.mat']);

sz=size(a.c);
nfound=0;
nchecked=0;
path=a.path;
pvs_acute=zeros(sz,'int16');  %the mask for PVS with acute angle pathes; 1: the acute angle path; 2: other PVS voxels
while nchecked<length(path)
    
    i=nchecked+1;
   
    pospath=path(i).subind(1:path(i).nvox_path,:);
    
    remove_path=false;
    for j=1:size(pospath,1)-2
        
        d1=pospath(j,:)-pospath(j+1,:);
        d2=pospath(j+1,:)-pospath(j+2,:);  %was j+2 and j+3 before; fixed 12/19/2018 
        cs=sum(d1.*d2)/sos(d1)/sos(d2);
        
        if cs<0
            pvs_acute(path(i).ind)=1;
            
            pvs_acute(pospath(j,1),pospath(j,2),pospath(j,3))=2;
            pvs_acute(pospath(j+1,1),pospath(j+1,2),pospath(j+1,3))=2;
            
            pvs_acute(pospath(j+2,1),pospath(j+2,2),pospath(j+2,3))=2;
            
            
          %{  
            nvox_path(i)=[];
            ind(i)=[];
            subind(i)=[];
            %}
            
            nfound=nfound+1;
            remove_path=true;
            [path1,path2]=break_pvs_acuteAngle(path(i),j+1);
            path(end+1:end+2)=[path1,path2];
            path(i)=[];
            break;
            
        end
        
    end
    %%
    if ~remove_path        
     nchecked=nchecked+1;
    end
end

if isfield(a,'history')

    history=a.history;
else
    history='';
end

if nfound>0
    mes=sprintf('%d paths found with >90 degree angle\n',nfound);
else
    mes='No paths with > 90 degree angle found.  Good!\n';
end

disp(mes);
history=[history,mes];

%orient=get_orient(prefix);
if nfound>0
 a.pvs_acute=pvs_acute;
end

[c,c_path]=set_c_from_ind(path,sz);
a.path=path;
a.history=history;
a.c=c;
a.cpath=c_path;
a.n_acuteAngle=nfound;
save(prefix,'-struct','a');


function [path1,path2]=break_pvs_acuteAngle(path,i_vertex) % break because of acute angle
% path contains ind, subind, and nvox_path
% i_vertex contains the index for the vertex of the acute angle

nvox_path=path.nvox_path;

xyz_path1=path.subind(1:i_vertex-1,:);
xyz_path2=path.subind(i_vertex+1:nvox_path,:);

ind_np=[path.ind(nvox_path+1:end),path.ind(i_vertex)]; %nonpath voxels
subind_np=[path.subind(nvox_path+1:end,:);path.subind(i_vertex,:)]; %nonpath voxels


%d1=distance_perpto_path(xyz_path1,xyz);
%d2=distance_perpto_path(xyz_path2,xyz);

[d1,gap1]=distance_to_path(xyz_path1,subind_np);
[d2,gap2]=distance_to_path(xyz_path2,subind_np);


path1.nvox_path=i_vertex-1;
path2.nvox_path=nvox_path-i_vertex;
sel2=(d1>=d2 & ~gap2)|(gap1&~gap2);
sel1=(d2>=d1 & ~gap1)|(gap2&~gap1);


path1.ind=[path.ind(1:i_vertex-1),ind_np(sel1)];
path2.ind=[path.ind(i_vertex+1:nvox_path),ind_np(sel2)];



path1.subind=cat(1,path.subind(1:i_vertex-1,:),subind_np(sel1,:));
path2.subind=cat(1,path.subind(i_vertex+1:nvox_path,:),subind_np(sel2,:));




 function res=distance_perpto_path(xyz_path,xyz)
  % distance of xyz to the line formed by principle component of the path
  % subind_path:l*3
  % subind: n*3
  % res: n*1
  
  k=pca(xyz_path);
  
  c=mean(xyz_path,1);
  res=0*xyz(:,1);
 
  for j=1:size(xyz,1)
     perp=xyz(j,:)-c-(xyz(j,:)-c)*k(:,1);    
     res(j)=sos(perp,2);     
  end


 function [res,gap]=distance_to_path(xyz_path,xyz)
  % minimum distance of xyz to path voxels xyz_path
     
  % subind_path:l*3
  % subind: n*3
  % res: n*1
  res=0*xyz(:,1);
  gap=res;
  
      for j=1:size(xyz,1)
          dist=sos(xyz(j,:)-xyz_path,2);
          [res(j),ind]=min(dist);
          pos=cat(1,xyz_path,xyz);
          gap(j)=hasgap(pos,ind,j+size(xyz_path,1));        
      end
  