dchi=linspace(0.2,1,100)*1e-6; 

m=0.44e-6*0.1*0.1;

for i=1:length(dchi)
s(i)=TotalSignal_Tissue_PerpPlane(sqrt(m/dchi(i)),dchi(i),7,10,90,1,0.3);

end
 
 figure;plot(dchi,s,'ro');