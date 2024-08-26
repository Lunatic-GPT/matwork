function aflip=flip_vector_origfirst(a)

aflip=zeros(size(a));
aflip(1)=a(1);

aflip(2:end)=a(end:-1:2);


%a=rand(1,128)+1i*rand(1,128);
%aflip=a(1);
%aflip(2:128)=a(128:-1:2);

%b=fft(aflip);

%b2=fft(a);
%b2flip=b2(1);
%b2flip(2:128)=b2(128:-1:2);
%b2flip is equal to aflip