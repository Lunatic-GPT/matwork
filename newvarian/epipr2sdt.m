
d2='Glu_Tau/Z_phantom_01';
dr='rawRE_noref';
di='rawIM_noref';

ns = readPar(d2,'ns');
nread=readPar(d2,'nread');
nphase=readPar(d2,'nphase');
nread=nread/2;
nnav=readPar(d2,'nnav');

str=dir(fullfile(dr,'*.fdf'));
if mod(length(str),ns)~=0
    error('number of files wrong');
end


ad=length(str)/ns;

br2=zeros(nread,nphase+nnav,ns,ad);
bi2=zeros(nread,nphase+nnav,ns,ad);

for i=1:ns
    for j=1:ad
    
     name = sprintf('slice%03dimage%03decho001.fdf',i,j);

     tmp=read_fdf(fullfile(dr,name));
     br2(:,:,i,j)=squeeze(tmp);
     
     name = sprintf('slice%03dimage%03decho001.fdf',i,j);

     tmp=read_fdf(fullfile(di,name));
     bi2(:,:,i,j)=squeeze(tmp);
     
    end
end

orient = readPar(d2,'orient');
orient=orient(2:end-1);

%if length(pss)>1 
 %   [pss_sort,isort]=sort(pss);
 %   b=b(:,:,isort,:);
%end

if strcmp(orient,'trans90')
    
  % b = permute(b,[2,1,3,4]);   
  b=flipdim(b,2);  %up down
   b=flipdim(b,1); % left right  %left on the left side in stimulate.
  ph=flipdim(ph,2);  %up down
   ph=flipdim(ph,1); % left right  %left on the left side in stimulate.
   
elseif strcmp(orient,'trans')
    b = permute(b,[2,1,3,4]); 
      b=flipdim(b,1);
      ph = permute(ph,[2,1,3,4]); 
      ph=flipdim(ph,1);
    
  %    b=flipdim(b,2);    
else
    b = permute(b,[2,1,3,4]);
      b=flipdim(b,1);
      b=flipdim(b,2);    %assume same as above, not tested.
    warning('The directions may be wrong');
end

if exist(d3,'dir')  
  writesdt4(ph,[d,'_iphs']);
end

if exist(d4,'dir')  
  writesdt4(br,[d,'_re']);
end
if exist(d5,'dir')  
  writesdt4(bi,[d,'_im']);
  writesdt4(angle(br+1i*bi),[d,'_angle']);
end


writesdt4(b,d);


