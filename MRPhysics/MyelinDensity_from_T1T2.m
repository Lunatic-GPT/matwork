R1s=1/1.07;
R1f=1/0.465;
R2s=1/0.117;
R2f=1/0.026;
f=[0.3,0.31];
TE = linspace(0,0.4);
TI= linspace(0,2);

[TE,TI]=meshgrid(TE,TI);

%% slow exchange limit
index=zeros([size(TE),length(f)]);

eindex=zeros([size(TE),length(f)]);

for i=1:length(f)
    S2 = f(i)*exp(-TE*R2f)+(1-f(i))*exp(-TE*R2s);
    S1 = f(i)*(1-2*exp(-TI*R1f))+(1-f(i))*(1-2*exp(-TI*R1s));
    index(:,:,i)=S1./S2;
    eindex(:,:,i)=sqrt(1./S1.^2+1./S2.^2);
end

%%

contrast=index(:,:,2)-index(:,:,1);
econtrast=sos(index.*eindex,3);

cnr=contrast./econtrast;

figure;imshow(cnr,[0,0.005]);


%% one sequence
TE = linspace(0,0.4);
TI= linspace(0,2);

[TE,TI]=meshgrid(TE,TI);
index=zeros([size(TE),length(f)]);

eindex=zeros([size(TE),length(f)]);

for i=1:length(f)
    S = f(i)*exp(-TE*R2f).*(1-2*exp(-TI*R1f))+(1-f(i))*exp(-TE*R2s).*(1-2*exp(-TI*R1s));
    
    index(:,:,i)=S;
    eindex(:,:,i)=1./abs(S);
end

contrast=index(:,:,2)-index(:,:,1);
econtrast=sos(index.*eindex,3);

cnr=contrast./econtrast;

figure;imshow(cnr*sqrt(2),[0,0.005]);



