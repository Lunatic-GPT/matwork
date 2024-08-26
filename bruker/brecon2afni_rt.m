function brecon2afni_rt(d,irecon)

if ~exist('irecon','var')
    irecon=1;
end

im=read_2dseq(d,irecon);
disp(size(im));
nr=readbPar(fullfile(d,'acqp'),'NR',true);
ni=readbPar(fullfile(d,'acqp'),'NI',true);
ns=readbPar(fullfile(d,'acqp'),'NSLICES',true);
im=reshape(im,[size(im,1),size(im,2),ns,length(im(:))/size(im,1)/size(im,2)/ns]);


if irecon==1
write_afni(im,d);
else
    write_afni(im,[d,'_',num2str(irecon)]);
end


function b=read_2dseq(scand,irecon)

if exist('2dseq','file')
    d=pwd;
elseif ~exist('irecon','var')
   d=fullfile(scand,'pdata','1');
else
    if ~isa(irecon,'char')
        irecon=num2str(irecon);
    end
   d=fullfile(scand,'pdata',irecon); 
end

endian=readbPar(fullfile(d,'reco'),'RECO_byte_order',false);

if strcmp(endian,'littleEndian')
 fid=fopen(fullfile(d,'2dseq'),'r','ieee-le');
else
 fid=fopen(fullfile(d,'2dseq'),'r','ieee-be');
end

sz=readbPar(fullfile(d,'reco'),'RECO_size');
tmp=readbPar(fullfile(d,'reco'),'RECO_offset');
transp=readbPar(fullfile(d,'reco'),'RECO_transposition');
ftsize=readbPar(fullfile(d,'reco'),'RECO_ft_size');

sz(3)=size(tmp,2);

%sz(1)=readbPar(fullfile(d,'d3proc'),'IM_SIX');
%sz(2)=readbPar(fullfile(d,'d3proc'),'IM_SIY');
%sz(3)=readbPar(fullfile(d,'d3proc'),'IM_SIZ');

fmt=readbPar(fullfile(d,'reco'),'RECO_wordtype',false);
switch fmt
    case '_32BIT_SGN_INT'
        
  %   d=fread(fid,sz(1)*sz(2)*sz(3),'int32');
     d=fread(fid,'int32');
  
    case '_16BIT_SGN_INT'
     %d=fread(fid,sz(1)*sz(2)*sz(3),'int16');
     
     d=fread(fid,'int16');
    case '_8BIT_UNSGN_INT'
     %  d=fread(fid,sz(1)*sz(2)*sz(3),'uint8');   
        d=fread(fid,'uint8');   
    case '_32BIT_FLOAT'
      %  d=fread(fid,sz(1)*sz(2)*sz(3),'float');
          d=fread(fid,'float');
    otherwise
end

 sz(4)=floor(length(d)/prod(sz));
 b=reshape(d(1:prod(sz)),sz');

if length(ftsize)==2
  for i=1:length(transp)
    if transp(i)==1
      tmp=reshape(b(:,:,i,:),[sz(2),sz(1),1,sz(4)]);
      b(:,:,i,:)=permute(tmp,[2,1,3,4]); 
    end
  end
else

 if transp(1)==1
  b=reshape(b,[sz(2),sz(1),sz(3),sz(4)]);
  b=permute(b,[2,1,3,4]);
 end
 
 
end
fclose(fid);
