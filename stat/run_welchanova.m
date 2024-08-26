function [fw,p]=run_welchanova(m)
% m: 2 D matrix with different columns corresponding to different groups

g=1:size(m,2);
g=repmat(g,[size(m,1),1]);

x=[m(:),g(:)];

[fw,p]=welchanova(x);
