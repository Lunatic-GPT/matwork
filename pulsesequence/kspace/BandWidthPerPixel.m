function bw=BandWidthPerPixel(bw,np,os)
% bw in Hz
% np: baseResolution;
% os: oversampling factor

dt=round(1e7/np/bw/os);  % in units of 100 ns

if dt>100
    dt = round(dt/100)*100;
else
    
    f100=[1,2,4,5,10,20,25,50,100]; % factors of 100
    
    for i=1:length(f100)-1
        
        if dt<=(f100(i)+f100(i+1))/2
            dt= f100(i);
            break;
        end
        
    end
    
    if dt>(f100(8)+f100(9))/2
        dt = 100;
    end
    
end

bw=1e7/np/dt/os;


if nargout==0
    disp(bw);
    disp(dt);
end