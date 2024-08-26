function dcm2pvs(dpattern,acpc,server,step_start)
% the acpc scan should be centered at the top of the corpus callosam.
%dcm2pvs(dpattern,acpc[,server,step_start])
% 10/21/2020: save ang.  Only calculate one lr position(i.e. one wip value).

if ~exist('server','var') || isempty(server)
    server='andrew';
end

if ~exist('dpattern','var')
    dpattern='tse_vfl_pss_FNRecon_matchR21*';
end


if ~exist('step_start','var')
    step_start=1;
end

th=tic;
dall=name4pat(dpattern);
dall=str2cell(dall);
dname=dall{1};

if step_start==1
    if ~exist([dname,'.nii.gz'],'file')
        dcm2nii(dname);
    end
end


if ~exist('acpc','var')
    acpc=name4pat('Slice_AC-PC*');
    acpc=str2cell(acpc);
    
    if ~isempty(acpc)
        acpc=acpc{1};
    end
    
end

if ~isempty(acpc)
    ori=readdPar(acpc,'ImageOrientationPatient');
    s_acpc.vec_lr=ori(1:3);
    norm=cross(ori(1:3),ori(4:6));
    [~,s_acpc.center]=dcmDimCenter(acpc);
    
    if norm(3)<0
        norm=-norm;
    end
    s_acpc.norm=norm;
    
end

if step_start<=2
    pvsseg_remote(dname,server);
end

if step_start==3
    pvsseg_remote(dname,server,1); %download only
end

if step_start<=4
    prob_prefix= pvsseg_remote(dname,server,2);
    prefix=gen_pvsmask(sprintf('%s.nii.gz',prob_prefix));
    process_pvs([prefix,'.nii.gz']);
else
    prefix=['mask_m2edn_nmsk_iter0_',dname];
end

if isempty(acpc)
    return;
end

load(fullfile(prefix,[prefix,'_stat.mat']),'pvs');

orient=get_orient_from_nii([dname,'.nii.gz']);

wip=3;
sgn=[-1,1,0];
label='rlb';

center=s_acpc.center+s_acpc.norm*15+s_acpc.vec_lr*20*sgn(wip);

res_search=Slice_P2PVS(pvs,orient,center,s_acpc.norm);

if ~isempty(res_search.th)
    pvs_perp=find_orient_max_PVS(res_search,10); %find the slice tilt with the most perpendicular PVS
end

save(sprintf('PCSlice_%s.mat',filename(dname)),'-struct','pvs_perp');
save(sprintf('PCSlice_%s.mat',filename(dname)),'-append','center','res_search','s_acpc');

fprintf('Perpendicular/Total PVSs: %d/%d \n',pvs_perp.n_perp,pvs_perp.n);

fprintf('%s - intersected pvs:',label(wip));
fprintf('perp - %d; total - %d ',pvs_perp.ipvs_perp);
fprintf('\nAngles: ');
fprintf('%3.1f ',pvs_perp.angles_perp);
fprintf('\n');


if ~exist('SlicePosition','dir')
    mkdir('SlicePosition');
end

fname=fullfile('SlicePosition',sprintf('SlicePosition%d_base.txt',wip));

save_slicePosition(fname,center,pvs_perp.norm,0);

fname2=fullfile('SlicePosition',sprintf('SlicePosition%d.txt',wip));

copyfile(fname,fname2);

toc(th);


function out=find_orient_max_PVS(res,thr)
% return the index array of intersecting PVSs and their tilt angles wrt the imaging plane

for i=1:length(res.angle)
    nperp(i)=sum(res.angle{i}<thr);
    n(i)=length(res.angle{i});
end
[nperp_max,imax]=max(nperp);
norm = res.norm(:,imax);
ipvs_perp=res.ipvs{imax}(res.angle{imax}<thr);

out.n_perp=nperp_max; % maximum number of perpendicular PVS
out.n = n(imax);
out.angles_perp=res.angle{imax}(res.angle{imax}<thr);
out.angles=res.angle{imax};
out.imax=imax;
out.thr=thr;
out.ipvs_perp=ipvs_perp;
out.ipvs=res.ipvs{imax};
out.norm=norm;


%%




