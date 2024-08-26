function  calc_DCF_truncate_proj( trj,thr,fname )
%calc_DCF_truncate_proj( trj,Nnew,fname )
% The second dimension is assumed the projection direction.
% Nnew is the new grid size

N=size(trj,2)*2;

trj2=sum(trj.^2,1);
trj2=sqrt(trj2);
sz=size(trj);

ind=length(find(trj2(1,:,1)<=thr/N/2));

trj=trj(:,1:ind,:,:)*N/thr;

DCF=calc_DCF(trj,thr);
DCF2=zeros([2,size(DCF)]);
DCF2(1,:,:,:,:)=DCF;
gdata=grid3_MAT_xp(DCF2,trj,ones([ind,sz(3:end)]),thr);
save(fname,'DCF','trj','gdata');

end
