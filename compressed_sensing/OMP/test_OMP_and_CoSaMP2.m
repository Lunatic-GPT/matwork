%{
Examples with Orthogonal Matching Pursuit (OMP) 
and Compressive Sampling Matching Pursuit (CoSaMP)

Also verifies that the codes are working on your computer;
see OMP.m and CoSaMP.m

Stephen Becker, Aug 1 2011
  small update, Oct 7 2011 (add the "oracle" solution for comparison)
%}

% -- Decide if we want real or complex-valued data
%COMPLEX = true;
 COMPLEX = false;

% -- Measurement matrix "A"

%% Add some noise:

x=zeros(1024,1);

id=randperm(1024);
x(id(1:100))=rand(100,1);
xr=x;
x=x+rand(1024,1)*0.1;

A=@(x) fft1c(x,1);
At=@(x) ifft1c(x,1);

b   = A(x);

    opts = [];
    opts.slowMode = 1;
    opts.printEvery     = 25;
    opts.maxiter=300;
    [xk] = OMP( {A,At}, b, 100,[],opts);
figure;plot(abs(xr),abs(x),'o');

xb=At(b);
figure;plot(abs(x),abs(xb),'o');

    