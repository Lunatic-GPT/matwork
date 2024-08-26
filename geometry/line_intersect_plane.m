function point_ints=line_intersect_plane(p1,p2,point,norm)
% point_ints=intersect_with_plane(p1,p2,point,norm)
% p1: 1*3
% p2: 1*3,
% point: a point in the plane
% norm: the normal direction of the plane
% point_ints: the intersect point of the line between p1 and p2 and the
% plane.  If not intersection, returns empty.

res=dist2plane([p1(:)';p2(:)'],point,norm);


if sign(res(1).*sign(res(2)))>0
    point_ints=[];
else
    point_ints= p1+(p2-p1)*res(1)/(res(1)-res(2));
end

