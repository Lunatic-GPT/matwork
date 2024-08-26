function suv=bqml2suv(val,dos,dname)
% suv=bqml2suv(val,dos,dname)
% or suv=bqml2suv(val,dname);
% val in units of bqml
% dos in mCi
% suv in kg/ltr
% dname: dicom file folder name
if nargin==2
 dname=dos;
 a=readdPar(dname,'RadiopharmaceuticalInformationSequence');
 dos=a.Item_1.RadionuclideTotalDose/3.7e7;
 
end

units=readdPar(dname,'Units');  % read dicom header
if ~strncmp(units,'BQML',4)
    disp(['Units: ', units]);
    error('Wrong units');
end



dos=dos*3.7e7; % convert to Bq

w=readdPar(dname,'PatientWeight');
w=w*1000;

df=readdPar(dname,'DecayFactor');
%sf=readdPar(dname,'ScatterFractionFactor');
sf=0;
suv=val*w*df/dos/(100-sf)*100;


function res3=readdPar(dname,par,allfile)
%res=readdPar(dname,par,allfile)

d=dir(fullfile(dname,'*.IMA'));
if isempty(d)
d=dir(fullfile(dname,'*.dcm'));

end

if ~exist('allfile','var') || ~allfile
 in=dicominfo(fullfile(dname,d(1).name));

 if exist('par','var')
  res3=getfield(in,par);
 else
    
    disp(in);
 end

else
    
    res3={};
    for i=1:length(d)
    in=dicominfo(fullfile(dname,d(i).name));

     res3{i}=getfield(in,par);
    end 
    
    
    if isa(res3{1},'char')
       return; 
    end
    
    res3=cell2mat(res3);
end

    