function fin_ROI=mdetect_kmeans(inner,outer,image,dark_roi)

n=1;


mn=mean(image(inner>0));
sd=std(image(inner>0));

for i=1:size(outer,1)
    for j=1:size(outer,2)
        if outer(i,j)>0 
           if  image(i,j)<mn/100 
             continue;
           end
          x(n)=image(i,j);
          iarr(n)=i;
          jarr(n)=j;
          n=n+1;
        end
    end
end


c1=mean(image(inner>0));
c2=mean(image(outer>0 & inner==0));

idx=kmeans(x,2,'Start',[c1;c2]);

fin_ROI=zeros(size(image));
for n=1:length(x)
    
    if idx(n)==1 
        fin_ROI(iarr(n),jarr(n))=1;
        
    end
    
    if ~any(inner & fin_ROI)
        fin_ROI=outer &~fin_ROI;
    end
    
end


fin_ROI=clusterize2(fin_ROI,4);

        

