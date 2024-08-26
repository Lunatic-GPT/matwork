function [m1,m2]=matrix_k2a(theta,phi)
% the matrix to convert from coefficients of spherical harmonics to those
% of polynomials along the six directions.

m1=zeros(length(theta),3);
m2=zeros(length(theta),5);
    
for i=1:length(theta)
       
       m_tmp=sp_harm(theta(i),phi(i));
       m1(i,:)=m_tmp(1:3);     
       m2(i,:)=m_tmp(4:8);
end

