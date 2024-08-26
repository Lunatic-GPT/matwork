function [a2,b2,iab,a_amb,b_amb,a_no,b_no]=clusters_overlap(a,b,max_clust_dist)
% assuming the clusters in a and b have different values
%a2 and b2 removes clusters that do not overlap but are otherwise the same as a and b
%iab: n*2, give the cluster voxel values for the overlapping clusters
%a_no: clusters that do not overlap with any cluster in b
%b_no: clusters that do not overlap with any cluster in a
%a_amb: ambiguous clusters in a that overlap with the same cluster in b
%b_amb: ambiguous clusters in b that overlap with the same cluster in a

ab=a+b;

ab=clusterize2(ab);

a2=a*0;
b2=b*0;

a_amb=a*0; %ambiguous clusters in a that overlap with the same cluster in b
b_amb=b*0; %ambiguous clusters in b that overlap with the same cluster in a


if ~exist('max_clust_dist','var')
    max_clust_dist=Inf;
end
ia=[];
ib=[];
na_2=0;  % clusters in a that overlap with multiple cluster in b
nb_2=0;  % clusters in b that overlap with multiple cluster in a
n=0;     % overlapping after excluding the above cases.

for i=1:max(ab(:))
    
    if mod(i,10)==0
     % disp([i,max(ab(:)),n]);
    end
    
    ma=ab(:)==i&a(:)>0;
    mb=ab(:)==i&b(:)>0;
    if any(ma) && any(mb)
        
        if isinf(max_clust_dist) || sos(roiCOM(ab==i&a>0) - roiCOM(ab==i&b>0)) <=max_clust_dist
          %  tic;
            tmpa=a(ma);
            tmpb=b(mb);
            
            unqa=unique(tmpa(:));
            unqb=unique(tmpb(:));
            unqa(unqa==0)=[];
            unqb(unqb==0)=[];
        %    toc;
            if length(unqa)>1
               nb_2=nb_2+1;
               
               for j=1:length(unqa)
                a_amb(a==unqa(j))=unqa(j);  
               end
               
               
            end
            
            if length(unqb)>1
               na_2=na_2+1; 
               for j=1:length(unqb)
                b_amb(b==unqb(j))=unqb(j);  
               end
               
            end
            
            if length(unqb)==1 && length(unqa)==1
                n=n+1;
                if n==100
                   disp(''); 
                end
                ia(n)=unqa;
                ib(n)=unqb;
                
                a2(ma)=unqa;          
                b2(mb)=unqb;
            end
         %   toc;
        end
    end
    
end

a_no = a;
a_no(a2>0|a_amb>0)=0;

b_no = b;
b_no(b2>0|b_amb>0)=0;

iab=[ia',ib'];
fprintf('%d clusters in a overlap with >1 cluster in b; excluded\n',na_2);
fprintf('%d clusters in b overlap with >1 cluster in a; excluded\n',nb_2);

fprintf('Total cluster (a/b/overlap) = (%d/%d/%d)\n',length(unique(a(:)))-1,length(unique(b(:)))-1,n);




