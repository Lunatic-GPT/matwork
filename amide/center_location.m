function res=center_location(dname)

%{
pos=readdPar(dname,'ImagePositionPatient');

nsl=double(readdPar(dname,'NumberOfSlices'));

row=double(readdPar(dname,'Rows'));
col=double(readdPar(dname,'Columns'));
dim=dcmdim(dname);

res=pos(1:3)'+[row,col,-(nsl-2)]/2.*dim(1:3);
%}
pos=readdPar(dname,'ImagePositionPatient');

row=double(readdPar(dname,'Rows'));
col=double(readdPar(dname,'Columns'));
dim=dcmdim(dname);
dim(3)=readdPar(dname,'SpacingBetweenSlices');
mat=[row-1,col-1,17];

res=pos(1:3)'+mat/2.*dim(1:3);






