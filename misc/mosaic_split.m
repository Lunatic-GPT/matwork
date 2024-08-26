function d3=mosaic_split(a,reps)


if length(reps)==1
    reps=[reps,reps];
end

sz=size(a);
if length(sz)==2
    sz(3)=1;
    sz(4)=1;
end

d=reshape(a,[sz(1)/reps(1),reps(1),sz(2)/reps(2),reps(2),sz(3:end)]);

d=permute(d,[1,3,4,2,5,6]);
d=reshape(d,[sz(1)/reps(1),sz(2)/reps(2),reps(1)*reps(2),sz(3:end)]);


d2=d(:,:,:,1);

d2=reshape(d2,[size(d2,1)*size(d2,2),size(d2,3)]);

sel=any(~isnan(d2)&d2~=0,1);

d3=d(:,:,sel,:);
