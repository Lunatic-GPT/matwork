function Mip_calc_afni( fname,ns,max_min,dim )
%Mip_calc( fname,ns,max_min,dim )


[d,info]=BrikLoad(fname);

if dim==1
    p=[3,2,1];
elseif dim==2
    
    p=[1,3,2];
else
    
    p=[1,2,3];
end

d2=permute(d,p);
d3=d2;
for i=1:size(d2,3)

    i1=floor(i-ns/2);
    i2=floor(i+ns/2-1);
    
    if i1<1
        i1=1;
    end
    if i2>size(d2,3)
        i2=size(d2,3);
    end
       
    
    if strcmp(max_min,'max')
     d3(:,:,i)=max(d2(:,:,i1:i2),[],3);
    else
     d3(:,:,i)=min(d2(:,:,i1:i2),[],3);
    end
        
    
    
end

d4=d3;
%d4=permute(d3,p);

prefix=strtok(fname,'.');

prefix=strtok(prefix,'+');
WriteBrikEZ(d4,info,'',sprintf('%s_MIP%d_dim%d',prefix,ns,dim),'');

