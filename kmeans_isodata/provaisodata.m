clear all
close all
clc
load('dades.mat')
ON=15;  %Threshold number of elements for the elimination of a grouping. 
OC=10;  %Threshold distance for the union of clusters.
OS=7;  %Standard deviation threshold for the division of a cluster
k=4;   %Number (maximum) clusters.
L=2;   %Maximum number of clusters can be mixed in a single iteration
I=10;  %Maximum number of iterations allowed.
NO=1;  %Extra parameter that does not automatically respond to a request for cambial some parameters..
min=50; %Minimum distance a point must be in each center. If you want to remove any item not give a high value. 
%%%%%%%%%%%%%%%%%%%%%
%  Funcion ISODATA  %
%%%%%%%%%%%%%%%%%%%%%
[centro, Xcluster, Ycluster, A, clustering]=isodata(X, Y, k, L, I, ON, OC, OS, NO, min);
clc;
fprintf('Numero de agrupaciones: %d',A);

% Presentacion de resultados por pantalla.

% Creamos los colores.
colr=zeros(A,3);
for i=1:A
    colr(i,:)=rand(1,3);
end;

% Representamos la informacion.
figure;
hold on;
for i=1:A,
    n=find(clustering==i);
    p=plot(X(n), Y(n),'.');set(p,'Color',colr(i,:));title(A)
end;

%plot(centro(:,1), centro(:,2), 'g.');

clc;
fprintf('Numero de agrupaciones: %d',A);
% Borramos variables temporales.
clear n;clear i;clear p;clear colr;clear NO;