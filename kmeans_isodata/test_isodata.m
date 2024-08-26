clear all
close all
clc

ON=15;  %Threshold number of elements for the elimination of a grouping. 
OC=10;  %Threshold distance for the union of clusters.
OS=7;  %Standard deviation threshold for the division of a cluster
k=10;   %Number (maximum) clusters.
L=3;   %Maximum number of clusters can be mixed in a single iteration
I=100;  %Maximum number of iterations allowed.
NO=1;  %Extra parameter that does not automatically respond to a request for cambial some parameters..
min=200; %Minimum distance a point must be in each center. If you want to remove any item not give a high value. 
%%%%%%%%%%%%%%%%%%%%%
%  Funcion ISODATA  %
%%%%%%%%%%%%%%%%%%%%%
temp = load('07.mat','image');

X =zeros(512*512,3);

for i=1:512
    for j=1:512
        d = sub2ind([512,512],i,j);
        X(d,1) = i;
        X(d,2) = j;
        X(d,3) = temp.image(i,j)/3;
    end
end
        
[centro, Xcluster, A, clustering]=isodata_ND(X, k, L, I, ON, OC, OS, NO, min);

Y = zeros(512,512);
for i=1:512*512
    [x1,x2] = ind2sub([512,512],i);
    Y(x1,x2) = clustering(i);
end
imagesc(Y);