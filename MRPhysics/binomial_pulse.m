function binomial_pulse()

pasRfSamples=zeros(2,220);

for lI=1:20
    pasRfSamples(1,lI)=0.5;
    pasRfSamples(2,lI)=0.0;
end
for lI=1:80
    pasRfSamples(1,20+lI)=0.0;
    pasRfSamples(2,20+lI)=0.0;
end

for lI=1:20
    pasRfSamples(1,100+lI)=1.0;
    pasRfSamples(2,100+lI)=180;
end

for lI=1:80
    pasRfSamples(1,120+lI)=0;
    pasRfSamples(2,120+lI)=0;
end

for lI=1:20
    pasRfSamples(1,200+lI)=0.5;
    pasRfSamples(2,200+lI)=0;
end

fprintf('%2.1f %d ',pasRfSamples(:));
fprintf('\n');