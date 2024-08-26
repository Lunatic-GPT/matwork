function cmat= remove_censor_points(mat,censor)

if length(censor) ~= size(mat,1)
    error('length(censor) ~= size(mat,1)');
end

clist = censor>0;
mat(clist,:) = [];
cmat = mat;

        