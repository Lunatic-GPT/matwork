


[k2,Data]=get_k_data(Data,Line,Partition,dfile,Nc,fov,lMatrix);
nufft=NUFFT3D(k2,1,[0,0,0],round(imSize/16),1,1);
im=nufft'*Data;
  
  
toc(tloc);


%{
adjoin:

imSize  nblock       estimated total time            per channel   nufft_init time
1/8      16          882                                0.5            37 
1/4      16          1100                               1              37
1         4          5150                                   36             137
1/4      8           1065                              1.7             75
1/16        1           





%}

% read 12.7 G data took 1085 s
% read 1 G data took 10 s


