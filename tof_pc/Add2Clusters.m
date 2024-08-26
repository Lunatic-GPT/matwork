function a_new=Add2Clusters(a,b)

%% add clusters in b to overlapping clusters in a; assign the same labels.

[a2,b2]=clusters_overlap(a,b);


b2val=unique(b2(:));
b2val=b2val(b2val~=0);

a_new=a;
for i=1:length(b2val)
  
  aval=a(a2==b2val(i));
  a_new(b2==b2val(i))=aval(1);


end
