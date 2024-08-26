function d=read_1dseq(scand,irecon)

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

endian=readbPar(fullfile(d,'visu_pars'),'VisuCoreByteOrder',false);

if strcmp(endian,'littleEndian')
 fid=fopen(fullfile(d,'2dseq'),'r','ieee-le');
else
 fid=fopen(fullfile(d,'2dseq'),'r','ieee-be');
end

sz=readbPar(fullfile(d,'reco'),'RECO_size');
transp=readbPar(fullfile(d,'reco'),'RECO_transposition');
ftsize=readbPar(fullfile(d,'reco'),'RECO_ft_size');

s=readbPar(fullfile(d,'reco'),'RECO_map_slope');
bl=readbPar(fullfile(d,'reco'),'RECO_map_offset');


%sz(3)=size(tmp,2);

%sz(1)=readbPar(fullfile(d,'d3proc'),'IM_SIX');
%sz(2)=readbPar(fullfile(d,'d3proc'),'IM_SIY');
%sz(3)=readbPar(fullfile(d,'d3proc'),'IM_SIZ');

fmt=readbPar(fullfile(d,'visu_pars'),'VisuCoreWordType',false);
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


fclose(fid);


