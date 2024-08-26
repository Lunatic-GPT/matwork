function ge3d2sdt(d)


d2=[d,'.fid'];
d1=[d,'.img'];
nread=readPar(d2,'np');
nphase=readPar(d2,'nv');
nphase2=readPar(d2,'nv2');
nread=nread/2;


str=dir(fullfile(d1,'*.fdf'));


    
b=read_fdf(fullfile(d1,'img_slab001image001echo001.fdf'));
     


orient = readPar(d2,'orient');
orient=orient(2:end-1);


switch orient
    case 'trans90'
    
  % b = permute(b,[2,1,3,4]);   
  b=flipdim(b,2);  %up down
   b=flipdim(b,1); % left right  %left on the left side in stimulate.

    case 'trans'
    b = permute(b,[2,1,3,4]); 
      b=flipdim(b,1);
  %    b=flipdim(b,2);    
    case 'sag'
    b = permute(b,[2,1,3,4]);
      b=flipdim(b,1);
      b=flipdim(b,2);    %assume same as above, not tested.
    case 'cor90'
    
    b=permute(b,    
    warning('The directions may be wrong');
end


 writesdt4(b,d);



