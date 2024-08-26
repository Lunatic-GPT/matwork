function [h,p,stats]=run_games_howell(m)

% m: 2 D matrix with different columns corresponding to different groups

g=1:size(m,2);
g=repmat(g,[size(m,1),1]);

[h,p,stats]=games_howell(m(:),g(:));
