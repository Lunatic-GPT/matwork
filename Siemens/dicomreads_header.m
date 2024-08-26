function res_struct=dicomreads_header(fname)



in=dicominfo(fname);

fn=fieldnames(in);
res_struct=struct();

for i=1:length(fn)
    d=getfield(in,fn{i});
   
    if isa(d,'uint8')  && strcmp(fn{i}(1:7),'Private')
        
        if any(d(1:8)~=['SV10',4,3,2,1]')
         continue;
        end
        
        element_count = char2long(d(9:12));
        
        current=17;
        
        for j=1:element_count
            name=sprintf('%c',d(current:current+63));
            name=strtok(name,' ');
            name=strrep(name,'-','_');
            subelement_count=char2long(d(current+76:current+79));
            current=current+84;
            values=cell(1,subelement_count);
            
            for k=1:subelement_count
                 
                len(1)=char2long(d(current:current+3));
                len(3)=char2long(d(current+8:current+11));
                len(4)=char2long(d(current+12:current+15));
                if len(1)~=len(4)
                    disp(fn{i});
                    disp('Wrong Format');
                end
                current=current+16;
                values{k}=sprintf('%c',d(current:current+len(1)-1));
                current=current+(ceil(len(1)/4)*4);    
                
            end
            
           res_struct.(name)=values;
        end
        
        
        
        return;
    end
 
 
end

function res=char2long(d)
d=double(d);
 res=d(1)+d(2)*256+d(3)*256^2+d(4)*256^3;

function res=us_get(dataset,tag,default)

if isfield(dataset,tag)
    res=getfield(data,tag);
else
    res=default;
end
