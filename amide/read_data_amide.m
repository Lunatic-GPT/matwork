function res = read_data_amide(fname,dname)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if ~exist('dname','var')
    dname=pwd;
end

xml1=[dname,'/data-set_',fname,'_raw-data.xml'];
format=readp(xml1,'raw_format');
dim=readp(xml1,'dim');

fid=fopen(sprintf('%s/data-set_%s_raw-data.dat',dname,fname));

  if strcmp(format, 'float-32-le')
        d=fread(fid,'float','ieee-le');
  elseif strcmp(format,'ubyte-8-ne')
      d=fread(fid,'uint8');
  elseif strcmp(format,'sshort-16-le')
      d=fread(fid,'int16');
  end

fclose(fid);



res=reshape(d,dim);


function res=readp(fname,par)

d=xmlread(fname);
thisListItem = d.getElementsByTagName(par);
childNode = thisListItem.item(0).getFirstChild;

res=str2num(childNode.getNodeValue);

if isempty(res)
    res=childNode.getNodeValue;
end

