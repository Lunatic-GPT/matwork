function data=myfft(data)

%    unsigned long n, mmax, m, j, istep, i;
%    double wtemp, wr, wpr, wpi, wi, theta;
%    double tempr, tempi;
 
 %   // reverse-binary reindexing
 nn=length(data)/2;
 
 n=nn*2;
 %   n = nn<<1;
    j=1;
    
    for i=1:2:n-1
        
        if (j>i) 
           
           tmp=data(i);
           data(i)=data(j);
           data(j)=tmp;
           
           tmp=data(i+1);
           data(i+1)=data(j+1);
           data(j+1)=tmp;
          
        end
        
        m = nn;
        while (m>=2 && j>m)
            j = j-m;
            m =floor(m/2);
        end
        
        j =j+ m;
    end
 
  
    mmax=2;
    while (n>mmax)
        istep = mmax*2;
        theta = -(2*pi/mmax);
        wtemp = sin(0.5*theta);
        wpr = -2.0*wtemp*wtemp;
        wpi = sin(theta);
        wr = 1.0;
        wi = 0.0;
        for m=1: 2:mmax-1
            
            for i=m:istep:n
                
                j=i+mmax;
                tempr = wr*data(j) - wi*data(j+1);
                tempi = wr * data(j+1) + wi*data(j);
                
                data(j) = data(i) - tempr;
                data(j+1) = data(i+1) - tempi;
                data(i) =data(i) + tempr;
                data(i+1) =data(i+1)+ tempi;
            end
            wtemp=wr;
            wr =wr+ wr*wpr - wi*wpi;
            wi = wi+wi*wpr + wtemp*wpi;
        end
        mmax=istep;
    end
    
    
    
    
    