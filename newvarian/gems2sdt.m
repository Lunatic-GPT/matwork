function gems2sdt(d)


d2=[d,'.fid'];
d1=[d,'.img'];
ns = readPar(d2,'ns');
nread=readPar(d2,'np');
nphase=readPar(d2,'nv');
nread=nread/2;

str=dir(fullfile(d1,'*.fdf'));
if mod(length(str),ns)~=0
    error('number of files wrong');
end

ad=length(str)/ns;

b=zeros(nread,nphase,ns,ad);
for i=1:ns
    for j=1:ad
    
     name = sprintf('slice%03dimage%03decho001.fdf',i,j);

     tmp=read_fdf(fullfile(d1,name));
     b(:,:,i,j)=squeeze(tmp);
    end
end



orient = readPar(d2,'orient');
orient=orient(2:end-1);
pss=readPar(d2,'pss');

if length(pss)>1 
    [pss_sort,isort]=sort(pss);
    b=b(:,:,isort,:);
end

if strcmp(orient,'trans90')
    
  % b = permute(b,[2,1,3,4]);   
  b=flipdim(b,2);  %up down
   b=flipdim(b,1); % left right  %left on the left side in stimulate.

elseif strcmp(orient,'trans')
    b = permute(b,[2,1,3,4]); 
      b=flipdim(b,1);
  %    b=flipdim(b,2);    
else
    b = permute(b,[2,1,3,4]);
      b=flipdim(b,1);
      b=flipdim(b,2);    %assume same as above, not tested.
    warning('The directions may be wrong');
end


    
writesdt4(b,d);



