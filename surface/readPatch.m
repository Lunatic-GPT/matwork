function [nodes,coords,tind,triangles]=readPatch(fname)
% [nodes,coords,tind,triangles]=readPatch(fname)
% reads an ascii flat patch file
% nodes
% coords: coordinates of the nodes
% tind: index of the triangles
% triangles:  the three node indices that make up a triangle.
  fid=fopen(fname,'r');
  cmt = fgetl(fid);
  if(cmt(1) ~='#')
      disp(['Wrong file format for a flat patch in', fname]);
      return;
  end

  counts = fscanf(fid,'%d',2);
  
  nodes = zeros(1,counts(1));
  coords = zeros(3,counts(1));
  tind = zeros(1,counts(2));
  triangles = zeros(3,counts(2));
  for nn = 1:counts(1)
      if feof(fid) 
          disp(['Wrong file format for a flat patch in', fname]);
          return;
      end
      
      ntemp =  fscanf(fid,'%d',1);
           
      fseek(fid,5,0);
      ntemp = fscanf(fid,'%d',1);
      nodes(nn) = ntemp;
      coords(:,nn) = fscanf(fid,'%f',3);
  end
  
  for nt = 1:counts(2)
      
      if feof(fid) 
          disp(['Wrong file format for a flat patch in', fname]);
          return;
      end
      tind(nt) = fscanf(fid,'%d',1);
      triangles(:,nt) = fscanf(fid,'%d',3);
  end
      